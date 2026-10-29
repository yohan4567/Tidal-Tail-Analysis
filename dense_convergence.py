"""Audit 8k-tracer real-input models without silently discarding spatial detail."""
import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from joint_grid import fixed_rate_fit,noise_report
from response_windows import rebin
from provenance import sha256_file


def thinning_multipliers(n,seed):
    """Nested unbiased one/two/four per arm per 1-Myr block from an 8k run.

    Complements/nested levels are not independent simulations.
    """
    if n<8 or n%8: raise ValueError('Two arms with complete four-particle blocks required')
    rng=np.random.default_rng(seed);rank=np.zeros(n,int)
    for start in range(0,n,4): rank[start:start+4]=rng.permutation(4)
    return {n//4:np.where(rank<1,4.,0.),n//2:np.where(rank<2,2.,0.),n:np.ones(n)}


def tv(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    if min(a.sum(),b.sum())<=0:return None
    return float(.5*np.abs(a/a.sum()-b/b.sum()).sum())


def analyse(study,selection,out):
    study,selection,out=map(Path,(study,selection,out))
    if out.exists():raise ValueError('Use a fresh convergence output')
    with np.load(selection/'maps.npz',allow_pickle=False) as d:
        counts=d['counts'][d['fitmask']];background=d['background'][d['fitmask']]
        x,y=np.nonzero(d['fitmask']); shape=d['fitmask'].shape;edges=d['edges']
        centre=(edges[:-1]+edges[1:])/2
        outer=(np.hypot(centre[x],centre[y])>=2)
    plan=json.loads((study/'plan.json').read_text())
    if plan['rates']!=[9.] or len(plan['seeds'])<3 or plan['release_interval_myr']!=.25:
        raise ValueError('This audit expects rate 9 and at least three 8k density realizations')
    factors=[1,2,4,8,16];times=sorted(plan['times_myr']);inputs={}
    seed_states={};normalization=[]
    for seed in plan['seeds']:
        run=study/f'rate_9_seed_{seed}'
        particles=pd.read_csv(run/'simulation/particles.csv')
        full=run/f'population_0_T_{max(times):g}'
        meta=json.loads((full/'response.json').read_text())
        path=full/'calibrated_counts.npz'
        if sha256_file(path)!=meta['template_sha256']:raise ValueError('Template hash mismatch')
        with np.load(path,allow_pickle=False) as d: fullstates=d['signal_draws']
        with np.load(full/'tracer_response.npz',allow_pickle=False) as d: tracer={k:d[k] for k in d.files}
        inputs[str(path.resolve())]=sha256_file(path);inputs[str((full/'tracer_response.npz').resolve())]=sha256_file(full/'tracer_response.npz')
        if len(particles)!=8000 or not np.array_equal(tracer['lookback_myr'],particles.lookback_release_myr.to_numpy()):
            # CSV decimal round trips can differ by a few ulps.
            if len(particles)!=8000 or not np.allclose(tracer['lookback_myr'],particles.lookback_release_myr,rtol=0,atol=1e-10):
                raise ValueError('8k particles and matched tracer ages required')
        rebuilt=rebin(tracer['calibration_counts'],tracer['fit_pixel_index'],tracer['lookback_myr'],max(times),len(counts))
        if not np.allclose(rebuilt,fullstates,rtol=1e-10,atol=1e-10):raise ValueError('Tracer response mismatch')
        levels=thinning_multipliers(len(particles),seed+10007)
        for density,multiplier in levels.items():
            for duration in times:
                states=rebin(tracer['calibration_counts']*multiplier[None,:],tracer['fit_pixel_index'],tracer['lookback_myr'],duration,len(counts))
                seed_states[(density,duration,seed)]=states
                mass=float(np.sum(tracer['mass_weight_msun']*multiplier*(tracer['lookback_myr']<=duration)))
                normalization.append(dict(seed=seed,density=density,duration_myr=duration,represented_mass=mass,passed=bool(np.isclose(mass,9*duration))))
    if not all(r['passed'] for r in normalization):raise ValueError('Thinned mass weights are not conserved')
    records=[];profiles=[];nested=[];pairs=[]
    for factor in factors:
        labels=(x//factor)*int(np.ceil(shape[1]/factor))+y//factor
        _,index=np.unique(labels,return_inverse=True);size=int(index.max()+1)
        def aggregate(a):
            if a.ndim==1:return np.bincount(index,weights=a,minlength=size)
            return np.array([np.bincount(index,weights=v,minlength=size) for v in a])
        n=aggregate(counts);b=aggregate(background)
        means={}
        for density in (2000,4000,8000):
            for duration in times:
                states=[aggregate(seed_states[(density,duration,seed)]) for seed in plan['seeds']]
                mean=np.mean(states,axis=0);means[(density,duration)]=mean
                beta,nll=fixed_rate_fit(n,mean,b)
                profiles.append(dict(factor=factor,density=density,T_myr=duration,seed='mean',nll=nll))
                if density==8000:
                    replicates=[s.mean(axis=0) for s in states]
                    noise=noise_report(replicates,mean.mean(axis=0)+beta*b)
                    outer_counts=np.array([seed_states[(density,duration,seed)].mean(axis=0)[outer].sum() for seed in plan['seeds']])
                    v=mean.mean(axis=0);order=np.argsort(-v);major=order[:np.searchsorted(np.cumsum(v[order]),.9*v.sum())+1]
                    se=np.std(replicates,axis=0,ddof=1)/np.sqrt(len(replicates))
                    lam=v+beta*b;ratios=se[major]/np.sqrt(lam[major])
                    records.append(dict(factor=factor,cell_deg=float((edges[1]-edges[0])*factor),T_myr=duration,bins=size,
                        expected_stars=float(v.sum()),major_pixel_target_passed=noise['major_pixel_target_passed'],
                        major_pixels_passing_fraction=float(np.mean(ratios<.1)),maximum_se_over_poisson=float(ratios.max()),
                        outer_region_se=float(outer_counts.std(ddof=1)/np.sqrt(len(outer_counts))),
                        outer_poisson_scale=float(np.sqrt(outer_counts.mean()+beta*background[outer].sum()))))
                    for i,j in combinations(range(len(states)),2):
                        pairs.append(dict(factor=factor,T_myr=duration,seed_a=plan['seeds'][i],seed_b=plan['seeds'][j],TV=tv(replicates[i],replicates[j])))
                    for seed,s in zip(plan['seeds'],states):
                        _,snll=fixed_rate_fit(n,s,b)
                        profiles.append(dict(factor=factor,density=8000,T_myr=duration,seed=str(seed),nll=snll))
        for lower,upper in ((2000,4000),(4000,8000)):
            for duration in times:
                a=means[(lower,duration)].mean(axis=0);v=means[(upper,duration)].mean(axis=0)
                nested.append(dict(factor=factor,T_myr=duration,lower_density=lower,upper_density=upper,TV=tv(a,v),
                    relative_count_change=float(abs(a.sum()-v.sum())/v.sum()),
                    scope='Nested thinning; correlated levels, not independent trajectories'))
    f=pd.DataFrame(records);p=pd.DataFrame(profiles)
    p['q']=2*(p.nll-p.groupby(['factor','density','seed']).nll.transform('min'))
    out.mkdir(parents=True)
    f.to_csv(out/'noise.csv',index=False);p.to_csv(out/'conditional_time_profiles.csv',index=False)
    pd.DataFrame(pairs).to_csv(out/'independent_seed_pairs.csv',index=False)
    pd.DataFrame(nested).to_csv(out/'nested_thinning.csv',index=False)
    result=dict(scope='Conditional rate 9/population 0 numerical test; no full MLR inference or physical validation.',
        normalization=normalization,noise=f.to_dict(orient='records'),
        fine_grid_noise_passed=int(f.loc[f.factor==1,'major_pixel_target_passed'].sum()),
        input_sha256=inputs,source_sha256=sha256_file(__file__),
        caveats=['Three seeds give uncertain variance estimates.',
            'Nested thinning is correlated and cannot independently prove convergence.',
            'Coarser grids lose spatial information; 1.78 deg cells are diagnostic only.',
            'Time profiles fix the physical rate and population, not profiles over those unknowns.'])
    (out/'convergence.json').write_text(json.dumps(result,indent=2))
    import os
    os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).parent/'outputs/.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for duration,g in f.groupby('T_myr'):
        axes[0].plot(g.cell_deg,g.maximum_se_over_poisson,'o-',label=f'{duration:g} Myr')
    axes[0].axhline(.1,color='black',ls='--');axes[0].set_xscale('log')
    axes[0].set(xlabel='Cell width (deg)',ylabel='Max SE / Poisson scale',title='8000 particles per seed, three seeds');axes[0].legend(fontsize=8)
    for density in (2000,4000,8000):
        g=p[(p.factor==1)&(p.seed=='mean')&(p.density==density)].sort_values('T_myr')
        axes[1].plot(g.T_myr/1000,g.q,'o-',label=f'{density} / seed')
    axes[1].set(xlabel='Lookback window (Gyr)',ylabel='Conditional q',title='Time profile at fixed rate and population');axes[1].legend(fontsize=8)
    fig.savefig(out/'dense_diagnostics.png',dpi=160);plt.close(fig)
    lines=['# 8000 粒子測試','',result['scope'],'',
        '| 像素寬度 deg | T Myr | MC 目標通過 | 最差 SE/Poisson |', '|---|---|---|---|']
    for row in records:lines.append(f"| {row['cell_deg']:.3f} | {row['T_myr']:g} | {row['major_pixel_target_passed']} | {row['maximum_se_over_poisson']:.3f} |")
    lines+=['','![Diagnostics](dense_diagnostics.png)','','詳細比較見 noise.csv、nested_thinning.csv、independent_seed_pairs.csv 與 conditional_time_profiles.csv。此測試不會自動升格成物理 MLR 測量。']
    (out/'REPORT_ZH.md').write_text('\n'.join(lines),encoding='utf-8')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('study','selection','out'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();analyse(a.study,a.selection,a.out)
