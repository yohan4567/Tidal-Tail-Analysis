"""Finite-grid Monte Carlo confidence sets, with matching production binning.

Coverage is conditional on frozen templates, calibration-state mixture and the
explicit background-scale grid. No claim between grid nodes or outside nuisances.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from joint_grid import read_templates,combined_states,aggregate_pixels,fixed_rate_fit
from profile_recovery import rank_threshold,statistics
from joint_recovery import wilson
from provenance import sha256_file


def profile_values(nll,keys):
    best=min(nll)
    rates=sorted({k[1] for k in keys}); times=sorted({k[0] for k in keys})
    return dict(rate_profile={r:max(0.,2*(min(v for v,k in zip(nll,keys) if k[1]==r)-best)) for r in rates},
        time_profile={t:max(0.,2*(min(v for v,k in zip(nll,keys) if k[0]==t)-best)) for t in times},
        time_rate_joint={(t,r):max(0.,2*(min(v for v,k in zip(nll,keys) if k[:2]==(t,r))-best))
                         for t in times for r in rates if any(k[:2]==(t,r) for k in keys)})


def project_set(values,accepted):
    values=sorted(values); accepted=sorted(accepted)
    return dict(accepted_grid_values=accepted,
        envelope=None if not accepted else [min(accepted),max(accepted)],
        touches_lower_grid_edge=bool(accepted and min(values) in accepted),
        touches_upper_grid_edge=bool(accepted and max(values) in accepted),
        scope='Discrete set only; envelope may contain rejected or unsampled values. No interpolated endpoints.')


def validate_grid(paths,scores_path,out,calibration=200,validation=300,seed=161803,background_scales=(1.,)):
    if calibration<20 or validation<20:raise ValueError('At least 20 trials per partition required')
    if not background_scales or any(not np.isfinite(b) or b<0 for b in background_scales):raise ValueError('Invalid background-scale grid')
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    n,b,groups,identity=read_templates(paths)
    score=json.loads(Path(scores_path).read_text(encoding='utf-8'))
    actual=sorted((str(Path(i['path']).resolve()),i['sha256']) for d in score['diagnostics'] for i in d['inputs'])
    expected=sorted((str(Path(p).resolve()),sha256_file(p)) for p in paths)
    if actual!=expected:raise ValueError('Scores and validation templates differ')
    labels=score.get('pixel_groups')
    keys=sorted(groups); signals=[combined_states(groups[k]) for k in keys]
    if labels is not None:
        labels=np.asarray(labels,int); n=aggregate_pixels(n,labels); b=aggregate_pixels(b,labels)
        signals=[aggregate_pixels(s,labels) for s in signals]
    observed=profile_values([fixed_rate_fit(n,s,b)[1] for s in signals],keys)
    for row in score['candidates']:
        if not np.isclose(observed['time_rate_joint'][(row['duration_myr'],row['rate_msun_myr'])],
                          min(r['q'] for r in score['candidates'] if (r['duration_myr'],r['rate_msun_myr'])==(row['duration_myr'],row['rate_msun_myr'])),atol=1e-6):
            raise ValueError('Saved grid likelihood does not match reconstructed binning')
    plan=dict(schema_version=1,calibration=calibration,validation=validation,seed=seed,
        background_scales=list(background_scales),score_sha256=sha256_file(scores_path),
        template_sha256=expected,source_sha256=sha256_file(__file__))
    plan_path=out/'plan.json'
    if plan_path.exists() and json.loads(plan_path.read_text())!=json.loads(json.dumps(plan)):
        raise ValueError('Resume plan differs; use a fresh directory')
    plan_path.write_text(json.dumps(plan,indent=2))
    cases=[]; thresholds={target:{} for target in observed}
    for i,truth in enumerate(keys):
        for j,beta in enumerate(background_scales):
            checkpoint=out/f'case_{i}_{j}.npz'
            if checkpoint.exists():
                with np.load(checkpoint,allow_pickle=False) as data: draws={k:data[k] for k in data.files}
            else:
                rng=np.random.default_rng(np.random.SeedSequence([seed,i,j])); draws={k:[] for k in observed}
                for trial in range(calibration+validation):
                    s=signals[i][rng.integers(len(signals[i]))]
                    counts=rng.poisson(s+beta*b)
                    stats=statistics([fixed_rate_fit(counts,v,b)[1] for v in signals],keys,truth)
                    for target in draws:draws[target].append(stats[target])
                draws={k:np.asarray(v) for k,v in draws.items()};np.savez(checkpoint,**draws)
            for target,values in draws.items():
                coordinate=truth[1] if target=='rate_profile' else truth[0] if target=='time_profile' else truth[:2]
                cutoff,rank=rank_threshold(values[:calibration])
                thresholds[target][coordinate]=max(thresholds[target].get(coordinate,0.),cutoff)
            cases.append((truth,beta,draws))
            print(f'Calibration truth {i+1}/{len(keys)}, background {beta:g} complete',flush=True)
    coverage=[]
    for truth,beta,draws in cases:
        row=dict(duration_myr=truth[0],rate_msun_myr=truth[1],scenario=truth[2],background_scale=beta)
        for target,values in draws.items():
            coordinate=truth[1] if target=='rate_profile' else truth[0] if target=='time_profile' else truth[:2]
            covered=values[calibration:]<=thresholds[target][coordinate]; count=int(covered.sum())
            row[target]=float(covered.mean());row[target+'_wilson95']=wilson(count,validation)
        coverage.append(row)
    accepted={target:[k for k,q in values.items() if q<=thresholds[target][k]] for target,values in observed.items()}
    report=dict(confidence_level=.95,method='Per-node Monte Carlo likelihood-ratio inversion; maximum cutoff over declared nuisance truths',
        mass_basis=identity[2],calibration_trials_per_truth=calibration,held_out_trials_per_truth=validation,
        coverage=coverage,rate=project_set(observed['rate_profile'],accepted['rate_profile']),
        time=project_set(observed['time_profile'],accepted['time_profile']),
        joint_accepted_grid_points=[list(k) for k in accepted['time_rate_joint']],
        validation_stage='pilot' if min(calibration,validation)<200 else 'conditional_finite_grid_calibration',
        sufficient_for_absolute_physical_rate=False,
        physical_rate_measurement_validated=False,
        limitations=['Coverage conditional on frozen finite template grid and explicit background-scale grid.',
            'Random calibration-state mixture; not coverage for every fixed calibration state.',
            'Finite Monte Carlo cutoffs vary; rank correction guarantees marginal calibration coverage, not conditional certainty for a realized bank.',
            'No between-node, off-grid, independent-template, mass-function or potential robustness guarantee.',
            'Particle refinement, outer-tail noise checks and background control-region validation remain required.'])
    if all(r['rate_profile']==1 for r in coverage):report['rate_uninformative_warning']='All held-out trials covered; grid may be uninformative or intervals conservative.'
    pd.DataFrame(coverage).to_csv(out/'coverage.csv',index=False)
    records=[]
    for (t,r),q in observed['time_rate_joint'].items():
        cutoff=thresholds['time_rate_joint'][(t,r)]
        records.append(dict(duration_myr=t,rate_msun_myr=r,q=q,cutoff95=None if np.isinf(cutoff) else cutoff,
                            cutoff_infinite=bool(np.isinf(cutoff)),accepted95=(t,r) in accepted['time_rate_joint']))
    pd.DataFrame(records).to_csv(out/'confidence_grid.csv',index=False)
    (out/'confidence95.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    from grid_plot import plot_grid
    plot_grid(score['candidates'],out/'confidence_heatmap.png',records)
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('templates-manifest','scores','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--calibration',type=int,default=200);p.add_argument('--validation',type=int,default=300)
    p.add_argument('--seed',type=int,default=161803)
    p.add_argument('--background-scales',type=float,nargs='+',default=[1.])
    a=p.parse_args();paths=json.loads(a.templates_manifest.read_text())['templates']
    validate_grid(paths,a.scores,a.out,a.calibration,a.validation,a.seed,a.background_scales)
