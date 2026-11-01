"""Weighted nuisance integration and relative likelihoods; no R/T posterior."""
import json,os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from provenance import sha256_file


def compare(grid,out,prior_path,beta_upper=5.):
    from bayesian_grid import _load_node,marginal_background
    grid,out,prior_path=map(Path,(grid,out,prior_path));score_path=grid/'scores/joint_grid.json'
    score=json.loads(score_path.read_text());prior=json.loads(prior_path.read_text());weights=prior['scenario_weights']
    import bayesian_grid,joint_fit
    plan=dict(method='nuisance_integrated_likelihood_no_target_posterior_v1',score_sha256=sha256_file(score_path),
        population_prior_sha256=sha256_file(prior_path),beta_weight=dict(distribution='uniform',lower=0.,upper=beta_upper),
        sources={str(Path(p).resolve()):sha256_file(p) for p in (Path(__file__),Path(bayesian_grid.__file__),Path(joint_fit.__file__))},
        rate_time_priors_used=False)
    if out.exists():
        if json.loads((out/'plan.json').read_text())!=plan:raise ValueError('Comparison plan changed; use a fresh output')
    else:out.mkdir(parents=True);(out/'plan.json').write_text(json.dumps(plan,indent=2))
    nodes=out/'integrated_nodes.jsonl';cached={}
    if nodes.exists():
        for line in nodes.read_text().splitlines():
            r=json.loads(line);cached[(r['scenario'],r['rate_msun_myr'],r['duration_myr'])]=r
    rows=[]
    for diagnostic in score['diagnostics']:
        n,s,b=_load_node(diagnostic,score.get('pixel_groups'))
        key=(diagnostic['scenario'],diagnostic['rate_msun_myr'],diagnostic['duration_myr'])
        if key in cached:r=cached[key]
        else:
            r=dict(scenario=key[0],rate_msun_myr=key[1],duration_myr=key[2],**marginal_background(n,s,b,beta_upper))
            with nodes.open('a',encoding='utf-8') as f:f.write(json.dumps(r)+'\n')
        rows.append(r)
    frame=pd.DataFrame(rows);result=[]
    for (rate,time),group in frame.groupby(['rate_msun_myr','duration_myr']):
        if set(group.scenario)!=set(weights):raise ValueError('Incomplete population quadrature at grid node')
        loglike=float(logsumexp([r.log_marginal_likelihood+np.log(weights[r.scenario]) for r in group.itertuples()]))
        result.append(dict(rate_msun_myr=float(rate),duration_myr=float(time),log_likelihood=loglike))
    f=pd.DataFrame(result);f['q']=2*(f.log_likelihood.max()-f.log_likelihood)
    f.to_csv(out/'likelihood_grid.csv',index=False)
    rp=f.groupby('rate_msun_myr',as_index=False).q.min();tp=f.groupby('duration_myr',as_index=False).q.min()
    rp.to_csv(out/'rate_profile.csv',index=False);tp.to_csv(out/'time_profile.csv',index=False)
    report=dict(plan=plan,best_sampled_node=f.loc[f.q.idxmin()].to_dict(),grid_points=len(f),
        target_parameter_posterior_computed=False,confidence_intervals_calibrated=False,physical_rate_measurement_validated=False,
        interpretation='Relative likelihood after background and measured-MF uncertainty integration; no prior or posterior over R/T. Delta scenarios remain separate.')
    (out/'comparison.json').write_text(json.dumps(report,indent=2))
    os.environ.setdefault('MPLCONFIGDIR',str(out/'mplcache'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from grid_plot import edges
    times=np.sort(f.duration_myr.unique());rates=np.sort(f.rate_msun_myr.unique());z=f.pivot(index='rate_msun_myr',columns='duration_myr',values='q').reindex(index=rates,columns=times).to_numpy()
    fig,ax=plt.subplots(1,3,figsize=(14,4),layout='constrained');im=ax[0].pcolormesh(edges(times)/1000,edges(rates),z,shading='flat',cmap='viridis_r');fig.colorbar(im,ax=ax[0],label='Relative -2 log likelihood')
    ax[0].set(xlabel='Lookback time (Gyr)',ylabel='Escape rate (solar masses/Myr)',title='Background and MF integrated')
    ax[1].plot(rp.rate_msun_myr,rp.q);ax[1].set(xlabel='Escape rate',ylabel='q',title='Profile over time')
    ax[2].plot(tp.duration_myr/1000,tp.q);ax[2].set(xlabel='Lookback time (Gyr)',ylabel='q',title='Profile over escape rate')
    fig.savefig(out/'likelihood_grid.png',dpi=160);plt.close(fig)
    return report
