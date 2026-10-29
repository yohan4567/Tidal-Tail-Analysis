"""Prespecified real-input joint-grid validation, with honest pass/fail reporting."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from joint_grid import read_templates,combined_states,analyse,noise_report,fixed_rate_fit
from joint_recovery import recover
from provenance import sha256_file
from validate_science import command,absolute_spec

HERE=Path(__file__).resolve().parent


def templates(folder):
    return [Path(p) for p in json.loads((folder/'templates.json').read_text())['templates']]


def mean_templates(paths):
    n,b,groups,identity=read_templates(paths)
    return n,b,{k:combined_states(v) for k,v in groups.items()},groups


def compare(a,b):
    sa=a.mean(axis=0); sb=b.mean(axis=0)
    return dict(expected_stars_a=float(sa.sum()),expected_stars_b=float(sb.sum()),
        relative_count_change=float(abs(sa.sum()-sb.sum())/max(sa.sum(),1e-30)),
        spatial_total_variation=float(.5*np.abs(sa/max(sa.sum(),1e-30)-sb/max(sb.sum(),1e-30)).sum()))


def run(out):
    out=Path(out).resolve()
    sources={k:sha256_file(HERE/k) for k in ('validate_joint_science.py','run_joint_study.py','joint_grid.py','joint_fit.py','joint_recovery.py',
        'response_model.py','response_windows.py','stream_model.py','dynamics.py','validate_window_reuse.py','mass_evolution.py')}
    plan=dict(source_sha256=sources,grid_rates=[4.,9.,25.],grid_times_myr=[250.,500.,750.,1000.],
        baseline_seeds=[42,84,126],baseline_interval_myr=5.,refined_seeds=[168,210,252],refined_interval_myr=2.5,
        short_seeds=[294,336,378],short_duration_myr=500.,short_interval_myr=2.5,
        populations=[-1,0,1],calibration_draws=2,observation_draws=16,recovery_trials=100,
        tolerances=dict(relative_count_change=.05,spatial_total_variation=.1,major_pixel_mc_fraction=.1),
        scope='Real Pal5 observations/calibration; conditional prescribed mass history; no independent N-body validation',
        physical_rate_measurement_validated=False)
    if out.exists():
        if not (out/'plan.json').exists() or json.loads((out/'plan.json').read_text())!=plan:
            raise ValueError('Validation plan/source differs; choose a fresh directory')
    else:
        out.mkdir(parents=True); (out/'plan.json').write_text(json.dumps(plan,indent=2))
    inputs=out/'inputs'; inputs.mkdir(exist_ok=True)
    config=json.loads((HERE/'outputs/evolving_mass_smoke_config.json').read_text())
    config['mass_evolution']['source_description']='Conditional Pal5 validation: f=1 living-star/dynamical escape conversion, fixed Plummer scale; no wind or remnant evolution'
    cfg=inputs/'dynamics.json'; cfg.write_text(json.dumps(config,indent=2))
    specs=[]
    for tag in ('minus1','plus0','plus1'):
        s=absolute_spec(HERE/f'outputs/Pal5_science_calibration/bundle/tilt_{tag}_zp_minus0.017.json')
        s.update(calibration_draws=2,observation_draws=16); s.pop('cache_dir',None)
        p=inputs/f'population_{tag}.json'; p.write_text(json.dumps(s,indent=2)); specs.append(p)
    selection=HERE/'outputs/Pal5_real_inputs_baseline/selection'
    # Freeze input hashes before simulations, including actual maps and calibrated inputs.
    input_files=[cfg,*specs,selection/'selection.json',selection/'maps.npz']
    for spec in specs:
        s=json.loads(spec.read_text())
        input_files += [Path(s['population'])/'population.csv',Path(s['population'])/'population.json',
            Path(s['quality'])/'quality.npz',Path(s['quality'])/'quality.json',Path(s['errors']['path']),Path(s['errors']['metadata']),Path(s['completeness']['path'])]
    hashes={str(p):sha256_file(p) for p in input_files}
    saved=out/'input_sha256.json'
    if saved.exists() and json.loads(saved.read_text())!=hashes: raise ValueError('Validation inputs changed')
    saved.write_text(json.dumps(hashes,indent=2))
    def study(name,rates,times,seeds,interval,populations):
        dest=out/name
        args=['--config',cfg,'--response-configs',*populations,'--selection',selection,'--out',dest,
            '--rates',*rates,'--times-myr',*times,'--seeds',*seeds,'--interval-myr',interval,'--cluster-age-myr',11000,'--workers',3]
        command('run_joint_study.py',args,out/(name+'.log'))
        print('Completed validation stage:',name,flush=True)
        return dest
    base=study('baseline',[4,9,25],[250,500,750,1000],[42,84,126],5,specs)
    refined=study('refined',[9],[250,500,750,1000],[168,210,252],2.5,[specs[1]])
    short=study('independent_short',[9],[500],[294,336,378],2.5,[specs[1]])
    # An off-grid physical rate receives its own dynamics, not an amplitude-scaled template.
    offgrid=study('off_grid_truth',[13],[375,1000],[420],5,[specs[1]])
    n,b,bmean,bgroups=mean_templates(templates(base))
    _,_,rmean,rgroups=mean_templates(templates(refined))
    _,_,smean,_=mean_templates(templates(short))
    comparisons=[]
    for key,states in rmean.items():
        row=compare(bmean[key],states); row.update(duration_myr=key[0],rate_msun_myr=key[1],scenario=key[2])
        row['passed_count_target']=row['relative_count_change']<=.05
        row['passed_shape_target']=row['spatial_total_variation']<=.1
        comparisons.append(row)
    shortkey=next(iter(smean)); shortcheck=compare(rmean[shortkey],smean[shortkey])
    shortcheck['passed_count_target']=shortcheck['relative_count_change']<=.05
    shortcheck['passed_shape_target']=shortcheck['spatial_total_variation']<=.1
    # Outer annulus predeclared by geometry, not by residuals.
    with np.load(selection/'maps.npz',allow_pickle=False) as d:
        centres=(d['edges'][1:]+d['edges'][:-1])/2
        x,y=np.meshgrid(centres,centres,indexing='ij'); outer=(np.hypot(x,y)>=2)[d['fitmask']]
    np.savez(out/'regions.npz',outer_annulus_ge_2deg=outer)
    if not (out/'region_scores').exists(): analyse(templates(base),out/'region_scores',out/'regions.npz')
    noise=[]
    for key,items in rgroups.items():
        states=rmean[key]; beta,_=fixed_rate_fit(n,states,b)
        noise.append(dict(duration_myr=key[0],report=noise_report([i['signal'].mean(axis=0) for i in items],states.mean(axis=0)+beta*b,{'outer_annulus_ge_2deg':outer})))
    # Use refined independent realizations as physical truth for on-grid rate 9 cases.
    def truth_path(folder,time,rate):
        return next(p for p in templates(folder) if f'_T_{time:g}' in str(p.parent.name) and f'rate_{rate:g}_' in str(p.parent.parent.name))
    scenarios=[
        dict(name='interior_time_sensitive',template=str(truth_path(refined,250,9)),tags=['interior','time_sensitive']),
        dict(name='interior_time_transition',template=str(truth_path(refined,500,9)),tags=['time_transition']),
        dict(name='long_window',template=str(truth_path(refined,1000,9)),tags=['boundary']),
        dict(name='low_signal',template=str(truth_path(base,250,4)),tags=['low_signal','boundary']),
        dict(name='background_only',null=True,tags=['null']),
        dict(name='physical_off_grid',template=str(truth_path(offgrid,375,13)),tags=['off_grid'])]
    manifest=out/'injections.json'
    manifest.write_text(json.dumps(dict(source_description='Prespecified frozen-real-input tests. Time sensitivity labels are provisional; saturation is assessed from measured contributions, not assumed.',scenarios=scenarios),indent=2))
    if not (out/'recovery/recovery.json').exists():
        recover(templates(base),manifest,out/'recovery',100)
    recovery=json.loads((out/'recovery/recovery.json').read_text())
    # Aggregate age contribution diagnostics before claiming saturation.
    age_tables=[]
    for p in templates(refined):
        if '_T_1000' in p.parent.name:
            d=pd.read_csv(p.parent/'age_contributions.csv'); d['source']=str(p.parent); age_tables.append(d)
    ages=pd.concat(age_tables,ignore_index=True); ages.to_csv(out/'age_preflight.csv',index=False)
    result=dict(plan=plan,refinement_comparisons=comparisons,independent_short_comparison=shortcheck,
        refined_noise=noise,recovery=recovery,
        physical_rate_measurement_validated=False,
        verdict='Completed conditional numerical/statistical validation experiment; inspect individual failures. Not a validation of total physical mass loss.',
        limits=['Finite 3-seed variance estimates; only one refined physical rate.',
            'No independent N-body injection or selection/background recalibration.',
            'No saturated-time scenario asserted without evidence; T tested only to 1 Gyr.',
            'Low-signal injection reuses a fitted template; other rate-9 truths use independent refined dynamics.'])
    (out/'validation.json').write_text(json.dumps(result,indent=2))
    lines=['# T–MLR validation with real Pal5 inputs','',
        'Completed experiment; individual checks may fail. No validated physical MLR is claimed.','',
        '## Refinement','', '| T (Myr) | Count change | Spatial TV | Count pass | Shape pass |','|---|---|---|---|---|']
    for row in comparisons:
        lines.append(f"| {row['duration_myr']:g} | {row['relative_count_change']:.2%} | {row['spatial_total_variation']:.3f} | {row['passed_count_target']} | {row['passed_shape_target']} |")
    lines+=['','## Independent short-window dynamics','',json.dumps(shortcheck,indent=2),'',
        '## Recovery','', 'Six scenarios, 100 trials each, split into calibration and held-out halves. See recovery/recovery.json for scenario-specific coverage, bias and missing scenarios.',
        '', '## Files','', '[Baseline grid](baseline/scores/REPORT.md) · [Region and noise diagnostics](region_scores/joint_grid.json) · [Full validation](validation.json)',
        '', 'Thresholds and sources were frozen in plan.json before computation. Failed convergence is reported, not hidden by changing tolerances.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    (out/'completion.json').write_text(json.dumps(dict(completed=True,source_sha256=sources,input_sha256=hashes),indent=2))
    print('Validation experiment complete:',out,flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--out',type=Path,required=True)
    run(p.parse_args().out)
