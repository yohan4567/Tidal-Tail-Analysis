"""Conditional injection/recovery on discrete joint templates, with held-out calibration.

This tests frozen templates/background, not selection or calibration construction.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from joint_grid import read_templates,fixed_rate_fit,combined_states
from provenance import sha256_file


def recover(paths,manifest,out,trials=100,seed=123):
    if trials<40: raise ValueError('At least 40 trials (calibration + held-out halves) required')
    out=Path(out)
    if out.exists(): raise ValueError('Use a fresh output directory')
    observed,b,groups,identity=read_templates(paths)
    candidates=[]
    for key,items in sorted(groups.items()):
        # Exactly the same seed-averaged templates used for the observed-data grid.
        candidates.append((key,combined_states(items)))
    spec=json.loads(Path(manifest).read_text(encoding='utf-8'))
    if not spec.get('source_description') or not spec.get('scenarios'):
        raise ValueError('Declared injection scenarios and source required')
    rng=np.random.default_rng(seed); records=[]; summaries=[]
    for case in spec['scenarios']:
        name=case['name']; beta=float(case.get('background_scale',1.))
        if beta<0 or not np.isfinite(beta): raise ValueError('Invalid injection background')
        truth_key=None; truth=None; truth_path=case.get('template')
        if case.get('null',False):
            truth=np.zeros((1,len(b)))
        else:
            if not truth_path: raise ValueError('Non-null injection requires independent/declared truth template')
            truth_path=(Path(manifest).resolve().parent/truth_path).resolve()
            _,tb,tg,ti=read_templates([truth_path])
            if ti!=identity or not np.array_equal(tb,b): raise ValueError('Injection observation/background mismatch')
            truth_key,items=next(iter(tg.items())); truth=items[0]['signal']
        ongrid=truth_key in groups if truth_key is not None else False
        case_records=[]
        for trial in range(trials):
            # Shared calibration state for a whole sky realization, not one draw per pixel.
            lam=truth[rng.integers(len(truth))]+beta*b
            counts=rng.poisson(lam)
            scores=[fixed_rate_fit(counts,s,b)[1] for _,s in candidates]
            best_index=int(np.argmin(scores)); best_key=candidates[best_index][0]
            _,null_nll=fixed_rate_fit(counts,np.zeros((1,len(b))),b)
            qnull=max(0.,2*(null_nll-min(scores)))
            qtruth=None
            if ongrid:
                # True parameters, evaluated using fitted template rather than injection realization.
                idx=[k for k,_ in candidates].index(truth_key)
                qtruth=max(0.,2*(scores[idx]-min(scores)))
            rec=dict(scenario=name,trial=trial,partition='calibration' if trial<trials//2 else 'held_out',
                fitted_duration_myr=best_key[0],fitted_rate_msun_myr=best_key[1],q_null=qnull,q_truth=qtruth)
            case_records.append(rec); records.append(rec)
        summary=dict(name=name,tags=case.get('tags',[]),null=bool(case.get('null',False)),on_grid=ongrid,
            truth_key=truth_key,truth_template_sha256=sha256_file(truth_path) if truth_path else None,
            conditional_scope='Frozen background/response templates; no selection recalibration.',trials=trials,
            coverage=None)
        cal=case_records[:trials//2]; test=case_records[trials//2:]
        if ongrid:
            cutoff=float(np.quantile([r['q_truth'] for r in cal],.95,method='higher'))
            covered=np.array([r['q_truth']<=cutoff for r in test])
            p=float(covered.mean())
            summary.update(empirical_95_q_threshold=cutoff,coverage=p,
                coverage_mc_standard_error=float(np.sqrt(p*(1-p)/len(test))),held_out_trials=len(test),
                rate_bias=float(np.mean([r['fitted_rate_msun_myr'] for r in test])-truth_key[1]))
            summary['coverage_wilson_95']=wilson(int(covered.sum()),len(covered))
        elif case.get('null',False):
            cutoff=float(np.quantile([r['q_null'] for r in cal],.95,method='higher'))
            exceed=np.array([r['q_null']>cutoff for r in test]); p=float(exceed.mean())
            summary.update(null_threshold=cutoff,held_out_false_positive_fraction=p,
                false_positive_mc_standard_error=float(np.sqrt(p*(1-p)/len(test))),held_out_trials=len(test))
            summary['false_positive_wilson_95']=wilson(int(exceed.sum()),len(exceed))
        else:
            summary.update(rate_bias=float(np.mean([r['fitted_rate_msun_myr'] for r in test])-truth_key[1]),
                coverage_caveat='Off-grid truth: no calibrated continuous confidence region is constructed.')
        summaries.append(summary)
    required={'interior','boundary','low_signal','null','time_sensitive','time_transition','time_saturated','off_grid'}
    supplied={tag for c in spec['scenarios'] for tag in c.get('tags',[])}
    result=dict(schema_version=1,seed=seed,trials_per_scenario=trials,scenarios=summaries,
        missing_required_scenario_tags=sorted(required-supplied),
        physical_measurement_validated=False,
        limitations=['Scenario tags are declarations, not automatic evidence of time sensitivity.',
            'Empirical thresholds are scenario-specific; not a universal observed-data confidence threshold.',
            'Finite frozen-template recovery does not validate population assumptions or independent N-body dynamics.'],
        input_sha256={str(Path(p).resolve()):sha256_file(p) for p in [*paths,manifest]})
    out.mkdir(parents=True); pd.DataFrame(records).to_csv(out/'trials.csv',index=False)
    (out/'recovery.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result


def wilson(successes,trials):
    z=1.959963984540054; p=successes/trials; d=1+z*z/trials
    centre=(p+z*z/(2*trials))/d
    half=z*np.sqrt(p*(1-p)/trials+z*z/(4*trials**2))/d
    return [max(0.,float(centre-half)),min(1.,float(centre+half))]


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--templates',type=Path,nargs='+',required=True); p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True); p.add_argument('--trials',type=int,default=100); p.add_argument('--seed',type=int,default=123)
    a=p.parse_args(); recover(a.templates,a.manifest,a.out,a.trials,a.seed)
