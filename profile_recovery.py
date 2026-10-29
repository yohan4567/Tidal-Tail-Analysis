"""Conditional recovery of profile confidence sets with finite calibration ranks.

Separates rate, time and time/rate statistics. Not a universal nuisance-robust belt.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from joint_grid import read_templates,combined_states,fixed_rate_fit
from joint_recovery import wilson
from provenance import sha256_file


def rank_threshold(scores,coverage=.95):
    values=np.asarray(scores,float)
    if values.ndim!=1 or not len(values) or not np.isfinite(values).all() or not 0<coverage<1:
        raise ValueError('Finite nonempty scores and coverage in (0,1) required')
    k=math.ceil((len(values)+1)*coverage)
    return (float(np.sort(values)[k-1]) if k<=len(values) else np.inf),k


def statistics(scores,keys,truth=None):
    scores=np.asarray(scores,float)
    if len(scores)!=len(keys) or not np.isfinite(scores).all(): raise ValueError('Finite grid scores required')
    best=float(scores.min()); result={}
    if truth is not None:
        masks=dict(rate_profile=[k[1]==truth[1] for k in keys],
            time_profile=[k[0]==truth[0] for k in keys],
            time_rate_joint=[k[:2]==truth[:2] for k in keys],
            candidate_point=[k==truth for k in keys])
        for name,mask in masks.items():
            if any(mask): result[name]=max(0.,2*(float(scores[np.asarray(mask)].min())-best))
    return result


def recover(paths,manifest,out,calibration=200,validation=300,seed=314159):
    if calibration<20 or validation<20: raise ValueError('At least 20 independent trials per partition required')
    out=Path(out)
    if out.exists(): raise ValueError('Use a fresh profile recovery directory')
    _,background,groups,identity=read_templates(paths)
    keys=sorted(groups); signals=[combined_states(groups[k]) for k in keys]
    config=json.loads(Path(manifest).read_text(encoding='utf-8'))
    if not config.get('source_description') or not config.get('scenarios'): raise ValueError('Declared scenarios required')
    out.mkdir(parents=True)
    rng=np.random.default_rng(seed); rows=[]; summaries=[]
    for scenario in config['scenarios']:
        name=scenario['name']; beta=float(scenario.get('background_scale',1.))
        if not np.isfinite(beta) or beta<0: raise ValueError('Invalid injected background scale')
        truth=None; truth_states=np.zeros((1,len(background)))
        if not scenario.get('null',False):
            if scenario.get('inject_exact_grid_mean',False):
                truth=(float(scenario['duration_myr']),float(scenario['rate_msun_myr']),json.dumps(scenario['scenario'],sort_keys=True))
                if truth not in groups: raise ValueError('Declared exact-grid truth missing')
                truth_states=signals[keys.index(truth)]
            else:
                path=(Path(manifest).resolve().parent/scenario['template']).resolve()
                _,b,tg,ti=read_templates([path])
                if ti!=identity or not np.array_equal(b,background): raise ValueError('Truth identity mismatch')
                truth,items=next(iter(tg.items()));truth_states=combined_states(items)
        case=[]
        for trial in range(calibration+validation):
            state=truth_states[rng.integers(len(truth_states))]
            counts=rng.poisson(state+beta*background)
            fitted=[fixed_rate_fit(counts,s,background)[1] for s in signals]
            _,null=fixed_rate_fit(counts,np.zeros((1,len(background))),background)
            scores=statistics(fitted,keys,truth)
            scores['null_search']=max(0.,2*(null-min(fitted)))
            chosen=keys[int(np.argmin(fitted))]
            row=dict(name=name,trial=trial,partition='calibration' if trial<calibration else 'held_out',
                fitted_rate=chosen[1],fitted_time=chosen[0],**scores)
            rows.append(row);case.append(row)
            if (trial+1)%50==0: print(f'{name}: {trial+1}/{calibration+validation}',flush=True)
        entry=dict(name=name,truth_key=truth,background_scale=beta,
            exact_grid_mean=scenario.get('inject_exact_grid_mean',False),statistics={})
        targets=['null_search'] if truth is None else [k for k in ('rate_profile','time_profile','time_rate_joint','candidate_point') if k in case[0]]
        for target in targets:
            threshold,k=rank_threshold([r[target] for r in case[:calibration]])
            held=np.array([r[target] for r in case[calibration:]])
            accepted=held<=threshold
            item=dict(threshold=float(threshold) if np.isfinite(threshold) else None,threshold_infinite=not np.isfinite(threshold),
                calibration_rank=k,calibration_trials=calibration,held_out_trials=validation,
                inclusion_fraction=float(accepted.mean()),inclusion_wilson_95=wilson(int(accepted.sum()),validation),
                expected_continuous_rank_coverage=k/(calibration+1) if k<=calibration else 1.)
            if target=='null_search':
                item.update(false_positive_fraction=float((~accepted).mean()),false_positive_wilson_95=wilson(int((~accepted).sum()),validation))
            entry['statistics'][target]=item
        summaries.append(entry)
        pd.DataFrame(rows).to_csv(out/'trials.csv',index=False)
        (out/'partial_summaries.json').write_text(json.dumps(summaries,indent=2))
    result=dict(schema_version=1,coverage_target=.95,seed=seed,scenarios=summaries,
        scope='Conditional, discrete grid, frozen background/selection and stated nuisance truths; not a uniform nuisance-robust confidence belt.',
        physical_rate_measurement_validated=False,
        input_sha256={str(Path(p).resolve()):sha256_file(p) for p in [*paths,manifest]},
        source_sha256=sha256_file(__file__))
    (out/'profile_recovery.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,required=True); p.add_argument('--templates-manifest',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--calibration',type=int,default=200);p.add_argument('--validation',type=int,default=300)
    p.add_argument('--seed',type=int,default=314159,help='Set an explicit different seed for independent follow-up experiments')
    a=p.parse_args();paths=json.loads(a.templates_manifest.read_text())['templates']
    recover(paths,a.manifest,a.out,a.calibration,a.validation,seed=a.seed)
