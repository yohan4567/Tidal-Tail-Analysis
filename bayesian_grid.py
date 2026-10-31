"""Grid posterior using raw Poisson mixtures and integrated background scale.

All intervals are Bayesian credible intervals conditional on explicit finite priors.
Never converts a profiled q curve directly into a marginalized posterior.
"""
import argparse,json,os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import logsumexp,xlogy
from joint_fit import fixed_rate_fit
from joint_grid import aggregate_pixels
from provenance import sha256_file


def integration_weights(x):
    x=np.asarray(x,float)
    if x.ndim!=1 or len(x)<2 or not np.isfinite(x).all() or (np.diff(x)<=0).any():raise ValueError('At least two strictly increasing finite axis nodes required')
    return .5*np.r_[x[1]-x[0],x[2:]-x[:-2],x[-1]-x[-2]]


def density_quantile(x,p,q):
    x=np.asarray(x,float);p=np.asarray(p,float)
    if p.shape!=x.shape or not 0<=q<=1 or not np.isfinite(p).all() or (p<0).any():raise ValueError('Invalid density or quantile')
    areas=.5*(p[:-1]+p[1:])*np.diff(x);total=areas.sum()
    if total<=0:raise ValueError('Empty density')
    if q==0:return float(x[0])
    if q==1:return float(x[-1])
    cumulative=np.r_[0.,np.cumsum(areas)];target=q*total;i=min(np.searchsorted(cumulative,target,side='right')-1,len(x)-2)
    width=x[i+1]-x[i];slope=(p[i+1]-p[i])/width
    u=brentq(lambda u:p[i]*u+.5*slope*u*u-(target-cumulative[i]),0,width)
    return float(x[i]+u)


def marginal_background(n,states,b,upper=5.):
    if not np.isfinite(upper) or upper<=0:raise ValueError('Positive prior beta upper bound required')
    n=np.asarray(n,float);states=np.atleast_2d(states);b=np.asarray(b,float)
    beta,nll=fixed_rate_fit(n,states,b)
    def scaled(v):
        lam=states+v*b
        ll=logsumexp(np.sum(xlogy(n,lam)-lam,axis=1))-np.log(len(states))
        return float(np.exp(ll+nll))
    width=max(.01,min(upper/4,5/np.sqrt(max(1.,n.sum()))))
    points=sorted({v for v in (beta-width,beta,beta+width) if 0<v<upper})
    def integrate(lo,hi,weighted=False):
        if hi<=lo:return 0.,0.
        return quad(lambda v:scaled(v)*(v if weighted else 1.),lo,hi,points=[v for v in points if lo<v<hi],
                    epsabs=1e-10,epsrel=1e-7,limit=200)
    area,error=integrate(0,upper)
    if not np.isfinite(area) or area<=0 or error>max(1e-10,area*1e-5):raise ValueError('Background integral unresolved; revise prior/support or numerical precision')
    moment,_=integrate(0,upper,True)
    lower,_=integrate(0,.05*upper);higher,_=integrate(.95*upper,upper)
    return dict(log_marginal_likelihood=float(-nll+np.log(area/upper)),background_mean=float(moment/area),
        beta_lower_edge_mass=float(lower/area),beta_upper_edge_mass=float(higher/area),quadrature_relative_error=float(error/area),profile_beta=beta)


def posterior_surface(nodes,population_weights,rate_prior='uniform',time_prior='uniform'):
    rates=np.array(sorted({r['rate_msun_myr'] for r in nodes}));times=np.array(sorted({r['duration_myr'] for r in nodes}))
    scenarios=sorted(population_weights);wr=integration_weights(rates);wt=integration_weights(times)
    if min(rates)<=0 or min(times)<=0:raise ValueError('Positive rate/time axes required')
    if rate_prior not in ('uniform','log_uniform') or time_prior not in ('uniform','log_uniform'):raise ValueError('Unknown axis prior')
    pr=np.ones(len(rates)) if rate_prior=='uniform' else 1/rates
    pt=np.ones(len(times)) if time_prior=='uniform' else 1/times
    pw=np.array([population_weights[s] for s in scenarios],float)
    if not np.isfinite(pw).all() or (pw<0).any() or pw.sum()<=0:raise ValueError('Invalid population prior')
    pw=pw/pw.sum();cube=np.full((len(scenarios),len(rates),len(times)),np.nan)
    ri={v:i for i,v in enumerate(rates)};ti={v:i for i,v in enumerate(times)};si={v:i for i,v in enumerate(scenarios)}
    for row in nodes:
        if row['scenario'] not in si:continue
        key=(si[row['scenario']],ri[row['rate_msun_myr']],ti[row['duration_myr']])
        if not np.isnan(cube[key]):raise ValueError('Duplicate posterior node')
        cube[key]=row['log_marginal_likelihood']
    if not np.isfinite(cube).all():raise ValueError('Population scenarios require the same complete rectangular rate/time grid')
    logdensity=logsumexp(cube+np.log(pw,where=pw>0,out=np.full_like(pw,-np.inf))[:,None,None],axis=0)+np.log(pr)[:,None]+np.log(pt)[None,:]
    shift=logdensity.max();density=np.exp(logdensity-shift);z=np.sum(density*wr[:,None]*wt[None,:]);density/=z
    rate_pdf=density@wt;time_pdf=wr@density
    def summarize(x,pdf):
        mass=pdf*integration_weights(x)
        return dict(median=density_quantile(x,pdf,.5),credible68=[density_quantile(x,pdf,.16),density_quantile(x,pdf,.84)],
            credible95=[density_quantile(x,pdf,.025),density_quantile(x,pdf,.975)],
            first_node_quadrature_mass=float(mass[0]),last_node_quadrature_mass=float(mass[-1]),
            sampled_mode_locations=[float(x[i]) for i in range(len(x)) if (i==0 or pdf[i]>=pdf[i-1]) and (i==len(x)-1 or pdf[i]>=pdf[i+1])],
            edge_diagnostic='Endpoint quadrature masses describe grid resolution and prior truncation; they are not tail probabilities beyond the grid.')
    summary=dict(rate_prior=rate_prior,time_prior=time_prior,rate=summarize(rates,rate_pdf),lookback_time=summarize(times,time_pdf),
        support=dict(rate_msun_myr=[float(rates.min()),float(rates.max())],time_myr=[float(times.min()),float(times.max())]),
        interpretation='Finite-support conditional Bayesian credible intervals, not frequentist confidence intervals or proof of physical accuracy.')
    return summary,dict(rates=rates,times=times,joint_density=density,rate_density=rate_pdf,time_density=time_pdf)


def _load_node(diag,pixel_groups):
    states=[];n=b=None
    for item in diag['inputs']:
        path=Path(item['path'])
        if sha256_file(path)!=item['sha256']:raise ValueError('Changed scored input: '+str(path))
        with np.load(path,allow_pickle=False) as d:
            nn=d['counts'].ravel();bb=d['background'].ravel();ss=d['signal_draws']
        if pixel_groups is not None:nn=aggregate_pixels(nn,pixel_groups);bb=aggregate_pixels(bb,pixel_groups);ss=aggregate_pixels(ss,pixel_groups)
        if n is None:n,b=nn,bb
        elif not np.array_equal(n,nn) or not np.array_equal(b,bb):raise ValueError('Mismatched observations/background')
        states.append(ss)
    return n,np.mean(states,axis=0),b


def run(grid,out,population_prior=None,beta_upper=5.):
    grid=Path(grid).resolve();out=Path(out).resolve();scorepath=grid/'scores/joint_grid.json';score=json.loads(scorepath.read_text())
    prior=json.loads(Path(population_prior).read_text()) if population_prior else None
    scenarios=sorted({d['scenario'] for d in score['diagnostics']})
    if prior and set(prior['scenario_weights'])!=set(scenarios):raise ValueError('Population-prior scenarios do not exactly match scored models')
    plan=dict(method='raw_poisson_mixture_background_integral_v1',scores_sha256=sha256_file(scorepath),
        source_sha256={p:sha256_file(Path(__file__).parent/p) for p in ('bayesian_grid.py','joint_fit.py','joint_grid.py')},
        population_prior_sha256=sha256_file(population_prior) if population_prior else None,
        beta_prior=dict(distribution='uniform',lower=0.,upper=beta_upper,chosen_without_fitting_to_rate=True),
        axis_priors=['uniform','log_uniform'],input_scope='Existing source-frozen grid only; no interpolation beyond support')
    if out.exists():
        if not (out/'plan.json').exists() or json.loads((out/'plan.json').read_text())!=plan:raise ValueError('Posterior plan differs; use a fresh directory')
    else:out.mkdir(parents=True);(out/'plan.json').write_text(json.dumps(plan,indent=2))
    if (out/'posterior_summary.json').exists():return json.loads((out/'posterior_summary.json').read_text())
    nodefile=out/'integrated_nodes.jsonl';cached={}
    if nodefile.exists():
        for line in nodefile.read_text().splitlines():
            r=json.loads(line);cached[(r['scenario'],r['rate_msun_myr'],r['duration_myr'])]=r
    nodes=[]
    for i,d in enumerate(score['diagnostics']):
        # Input hashes are checked even when using a cached background integral.
        key=(d['scenario'],d['rate_msun_myr'],d['duration_myr'])
        for item in d['inputs']:
            if sha256_file(item['path'])!=item['sha256']:raise ValueError('Changed scored input: '+item['path'])
        if key in cached:r=cached[key]
        else:
            n,s,b=_load_node(d,score.get('pixel_groups'));r=dict(scenario=d['scenario'],rate_msun_myr=d['rate_msun_myr'],duration_myr=d['duration_myr'],**marginal_background(n,s,b,beta_upper))
            with nodefile.open('a') as f:f.write(json.dumps(r)+'\n');f.flush()
        nodes.append(r)
        if (i+1)%50==0:print(f'Background integrals {i+1}/{len(score["diagnostics"])}',flush=True)
    results={};surfaces={}
    models={'measured_mf_marginalized':prior['scenario_weights']} if prior else {f'conditional_scenario_{i}':{s:1.} for i,s in enumerate(scenarios)}
    for model,weights in models.items():
        for rp,tp in [('uniform','uniform'),('log_uniform','uniform'),('uniform','log_uniform'),('log_uniform','log_uniform')]:
            name=f'{model}__rate_{rp}__time_{tp}';summary,arrays=posterior_surface(nodes,weights,rp,tp)
            summary['population_weights']=weights;results[name]=summary;surfaces[name]=arrays
            np.savez_compressed(out/(name+'.npz'),**arrays)
    pd.DataFrame(nodes).to_csv(out/'integrated_nodes.csv',index=False)
    os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).parent/'outputs/.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for name,a in surfaces.items():
        axes[0].plot(a['rates'],a['rate_density'],label=name.replace('conditional_scenario_','scenario '),alpha=.8)
        axes[1].plot(a['times']/1000,a['time_density']*1000,alpha=.8)
    axes[0].set(xlabel='Escape rate (solar masses/Myr)',ylabel='Conditional posterior density',title='Background and time integrated')
    axes[1].set(xlabel='Lookback time (Gyr)',ylabel='Conditional posterior density / Gyr',title='Background and rate integrated')
    axes[0].legend(fontsize=5);fig.savefig(out/'posterior_prior_sensitivity.png',dpi=150);plt.close(fig)
    primary=next(iter(surfaces));a=surfaces[primary]
    fig,ax=plt.subplots(figsize=(7,5),layout='constrained');im=ax.pcolormesh(a['times']/1000,a['rates'],a['joint_density'],shading='auto');fig.colorbar(im,ax=ax,label='Joint posterior density / (Msun/Myr * Myr)')
    ax.set(xlabel='Lookback time (Gyr)',ylabel='Escape rate (solar masses/Myr)',title='Conditional joint posterior; finite grid support');fig.savefig(out/'joint_posterior.png',dpi=150);plt.close(fig)
    report=dict(method=plan['method'],plan=plan,population_combined=prior is not None,
        population_prior=prior,results=results,physical_rate_measurement_validated=False,
        credible_intervals_are_confidence_intervals=False,finite_support=True,
        maximum_background_quadrature_relative_error=max(r['quadrature_relative_error'] for r in nodes),
        limitations=['All credible intervals are conditional on explicit finite rate/time and background priors.',
            'Uniform vs log-uniform sensitivity is reported; reaching a support edge does not establish a physical bound.',
            'No equal-probability mixing of arbitrary legacy tilt models: without a population prior, report each separately.',
            'New measured-MF priors constrain progenitor global MF; escape-transfer assumptions remain separate scenarios.',
            'Finite particle/calibration quadrature and population-node convergence require validation.',
            'Adopted orbit, halo, background shape, extinction and unseen stellar masses remain conditional.'])
    lines=['# Conditional Bayesian grid inference','','Background normalization is integrated from raw Poisson-mixture templates; profiled q is not treated as a posterior.','',
        'All intervals below are finite-support Bayesian credible intervals, not calibrated frequentist confidence intervals.','',
        '| Model / priors | MLR median | MLR 68% | MLR 95% | T median (Myr) |','|---|---:|---|---|---:|']
    for name,r in results.items():lines.append(f"| {name} | {r['rate']['median']:.3f} | {r['rate']['credible68']} | {r['rate']['credible95']} | {r['lookback_time']['median']:.1f} |")
    lines+=['','![Prior sensitivity](posterior_prior_sensitivity.png)','','![Joint posterior](joint_posterior.png)','',*report['limitations']]
    (out/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')
    (out/'posterior_summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--grid',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--population-prior',type=Path);p.add_argument('--beta-upper',type=float,default=5.);a=p.parse_args()
    run(a.grid,a.out,a.population_prior,a.beta_upper)
