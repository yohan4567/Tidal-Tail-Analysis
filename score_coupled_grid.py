"""Score rate-coupled dynamics at their physical rate: profile background only."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from fit_rate import ensemble_objective,profile_background
from provenance import sha256_file


def load_template(path):
    path=Path(path)
    with np.load(path,allow_pickle=False) as d:
        n=d['counts']; b=d['background']; ref=float(d['reference_mdot'])
        s=d['signal_draws'] if 'signal_draws' in d else d['signal'].ravel()[None,:]
        desc=json.loads(str(d['selection_description'].item()))
    identity={k:desc.get(k) for k in ('observation_identity','selection_sha256','mass_basis')}
    if not identity['observation_identity'] or not np.isfinite(ref) or ref<=0:
        raise ValueError('Structured observation identity and positive reference rate required')
    return n,s,b,ref,desc,identity


def score(paths,out,*,fixed_templates=()):
    if not paths: raise ValueError('At least one coupled candidate required')
    out=Path(out)
    if out.exists(): raise ValueError('Use a fresh grid score directory')
    reference=None; identity=None; rows=[]; rates=set()
    for path in map(Path,paths):
        n,s,b,ref,desc,key=load_template(path)
        history=desc.get('dynamics_mass_history')
        if desc.get('amplitude_policy')!='fixed_reference_rate' or history is None or not np.isclose(ref,history['stellar_rate_msun_myr']):
            raise ValueError('Coupled template reference rate and prescribed mass history must agree')
        if reference is None: reference=(n,b); identity=key
        elif key!=identity or not np.array_equal(n,reference[0]) or not np.array_equal(b,reference[1]):
            raise ValueError('All candidates need identical observations/background/mass basis')
        objective,_,_,_,_=ensemble_objective(n,s,b)
        background,nll=profile_background(objective,1.,n.sum()/b.sum())
        rates.add(ref)
        rows.append(dict(dynamics='evolving_mass',rate_msun_myr=ref,population_tilt=desc['scenario']['escape_mass_tilt'],
            background_scale=background,negative_log_likelihood=nll,signal_amplitude=1.,
            template=str(path.resolve()),template_sha256=sha256_file(path)))
    for path in map(Path,fixed_templates):
        n,s,b,ref,desc,key=load_template(path)
        if key!=identity or not np.array_equal(n,reference[0]) or not np.array_equal(b,reference[1]) or desc.get('dynamics_mass_history') is not None:
            raise ValueError('Fixed controls must match the same dataset and have rate-independent dynamics')
        objective,_,_,_,_=ensemble_objective(n,s,b)
        for rate in sorted(rates):
            background,nll=profile_background(objective,rate/ref,n.sum()/b.sum())
            rows.append(dict(dynamics='fixed_mass',rate_msun_myr=rate,population_tilt=desc['scenario']['escape_mass_tilt'],
                background_scale=background,negative_log_likelihood=nll,signal_amplitude=rate/ref,
                template=str(path.resolve()),template_sha256=sha256_file(path)))
    frame=pd.DataFrame(rows)
    frame['delta_minus2_log_likelihood']=2*(frame.negative_log_likelihood-frame.negative_log_likelihood.min())
    out.mkdir(parents=True)
    frame.to_csv(out/'scores.csv',index=False)
    result=dict(observation_identity=identity,candidates=frame.to_dict(orient='records'),
        interpretation='Discrete physically coupled rate candidates. Only background normalization is fitted for each evolving candidate; no free signal amplitude or interpolated confidence interval.',
        limitations=['Prescribed mass history with assumed dynamical/population mass conversion, not self-consistent internal evolution.',
            'Population/radius/Galactic potential/extinction assumptions remain conditional.',
            'Coarse grid and finite tracers do not determine a continuous best rate or calibrated model-selection significance.'])
    (out/'scores.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parent/'outputs/.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
    colours={-1:'#2462a8',0:'#d97922',1:'#2d9859'}
    for (kind,tilt),group in frame.groupby(['dynamics','population_tilt']):
        group=group.sort_values('rate_msun_myr')
        ax.plot(group.rate_msun_myr,group.delta_minus2_log_likelihood,
            'o-' if kind=='evolving_mass' else 's--',color=colours.get(tilt),
            label=f"{'Evolving' if kind=='evolving_mass' else 'Fixed'} mass / population tilt {tilt:+g}")
    ax.set(xlabel='Physical assumed escape rate (solar masses / Myr)',ylabel='Relative -2 log likelihood',
        title='Rate changes also change the dynamics in evolving-mass candidates')
    ax.legend(fontsize=8); ax.grid(alpha=.2)
    fig.suptitle('Lines connect sampled rates only; no continuous optimum or confidence interval',fontsize=10)
    fig.savefig(out/'coupled_grid.png',dpi=160); plt.close(fig)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--templates',type=Path,nargs='+',required=True)
    p.add_argument('--fixed-templates',type=Path,nargs='*',default=[]); p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); r=score(a.templates,a.out,fixed_templates=a.fixed_templates)
    print('Scored candidates:',len(r['candidates']))
