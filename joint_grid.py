"""Conditional T/rate grid, trajectory-window reuse and independent-template diagnostics.

No asymptotic confidence intervals or automatic physical-validation claim.
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from joint_fit import fixed_rate_fit

from provenance import sha256_file


def window_weights(particles, run, duration, weights):
    """Keep original weights; exact stratified normalization only on bin boundaries."""
    maximum=float(run['config']['duration_myr']); duration=float(duration)
    if not np.isfinite(duration) or not 0<duration<=maximum:
        raise ValueError('Window must be finite and within simulated lookback range')
    ages=particles.lookback_release_myr.to_numpy(float)
    w=np.asarray(weights,float).copy()
    if len(w)!=len(ages) or not np.isfinite(ages).all() or (ages<0).any() or (ages>maximum+1e-7).any():
        raise ValueError('Invalid release ages or weight shape')
    if duration<maximum:
        history=run.get('dynamics_mass_history')
        if not history or history.get('mode')!='prescribed_linear_escape':
            raise ValueError('Short windows currently require a prescribed constant-rate coupled history')
        sampling=run.get('release_time_sampling')
        if sampling not in ('uniform_lookback_both_arms','stratified_lookback_both_arms'):
            raise ValueError('Unsupported release sampling for subset reuse')
        if sampling=='stratified_lookback_both_arms':
            width=maximum/(run['config']['n_particles']/2)
            if not np.isclose(duration/width,round(duration/width),rtol=0,atol=1e-7):
                raise ValueError('Stratified window must end on a release-bin boundary')
    w[ages>duration]=0.
    return w,dict(duration_myr=duration,maximum_duration_myr=maximum,
        reuse='release_age_subset_original_weights',retained_tracers=int(np.count_nonzero(w)),
        represented_mass_msun=float(w.sum()),
        expected_mass_msun=float(run['config']['reference_mdot_msun_myr']*duration),
        caveat='Not a cluster formation time. Random uniform subsets have finite-sampling mass fluctuations; never renormalize them.')


def noise_report(replicates, expected, regions=None):
    """Independent trajectory realizations; calibration states are NOT replicates."""
    a=np.asarray(replicates,float); lam=np.asarray(expected,float)
    if a.ndim!=2 or a.shape[1:]!=lam.shape or len(a)<2:
        return dict(status='unverified',reason='At least two independent dynamics realizations required')
    se=a.std(axis=0,ddof=1)/np.sqrt(len(a)); mean=a.mean(axis=0)
    order=np.argsort(-mean); cumulative=np.cumsum(mean[order])
    major=order[:np.searchsorted(cumulative,.9*mean.sum())+1] if mean.sum()>0 else np.array([],int)
    passed=bool(len(major) and np.all(se[major]<.1*np.sqrt(lam[major])))
    result=dict(status='diagnostic_only',independent_realizations=len(a),major_pixel_target_passed=passed,
        target='SE(mean template) < 0.1 sqrt(expected total counts), top 90% signal-count pixels',
        maximum_standard_error=float(se.max()),regions={},
        caveat='Few seeds give uncertain variance. Also require density/refinement and likelihood stability; no automatic convergence certification.')
    for name,mask in (regions or {}).items():
        mask=np.asarray(mask,bool)
        if mask.shape!=lam.shape or not mask.any(): raise ValueError('Invalid predeclared region mask')
        values=a[:,mask].sum(axis=1) # preserves spatial covariance
        result['regions'][name]=dict(mean=float(values.mean()),standard_error=float(values.std(ddof=1)/np.sqrt(len(a))),
            poisson_scale=float(np.sqrt(lam[mask].sum())))
    return result


def read_templates(paths):
    groups={}; identity=None; observed=None; background=None
    for path in map(Path,paths):
        with np.load(path,allow_pickle=False) as d:
            desc=json.loads(str(d['selection_description'].item()))
            n=d['counts'].ravel(); b=d['background'].ravel()
            s=d['signal_draws'] if 'signal_draws' in d else d['signal'].ravel()[None,:]
            rate=float(d['reference_mdot'])
        if desc.get('amplitude_policy')!='fixed_reference_rate' or not desc.get('lookback_window'):
            raise ValueError('Versioned, windowed coupled response templates required')
        history=desc.get('dynamics_mass_history') or {}
        if not np.isclose(rate,history.get('stellar_rate_msun_myr',np.nan)):
            raise ValueError('Rate/history mismatch')
        key=tuple(desc.get(k) for k in ('observation_identity','selection_sha256','mass_basis','dynamics_family_identity','model_version'))
        if not all(key): raise ValueError('Missing observation identity/mass basis/model version')
        if identity is None: identity=key; observed=n; background=b
        if key!=identity or not np.array_equal(n,observed) or not np.array_equal(b,background):
            raise ValueError('Grid observations, background, mass basis and model version must match')
        scenario=json.dumps(desc.get('scenario'),sort_keys=True)
        group=(float(desc['lookback_window']['duration_myr']),rate,scenario)
        groups.setdefault(group,[]).append(dict(signal=s,path=str(path.resolve()),sha256=sha256_file(path),desc=desc))
    if not groups: raise ValueError('No templates supplied')
    return observed,background,groups,identity


def combined_states(items):
    shapes={i['signal'].shape for i in items}
    if len(shapes)!=1: raise ValueError('Independent templates need aligned calibration-state dimensions')
    response_ids=[i['desc'].get('response_sampling_identity') for i in items]
    dynamics_ids=[i['desc'].get('dynamics_sampling_identity') for i in items]
    if len(items)>1 and (None in response_ids or len(set(response_ids))!=1 or None in dynamics_ids or len(set(dynamics_ids))!=len(items)):
        raise ValueError('Independent dynamics realizations and aligned response sampling identities required')
    return np.mean([i['signal'] for i in items],axis=0)


def aggregate_pixels(values, index):
    """Sum disjoint pixels, preserving counts and calibration-state alignment."""
    values=np.asarray(values); index=np.asarray(index)
    if index.ndim!=1 or len(index)!=values.shape[-1] or not np.issubdtype(index.dtype,np.integer) or (index<0).any():
        raise ValueError('Nonnegative integer labels required for every fitted pixel')
    _,index=np.unique(index,return_inverse=True)
    flat=values.reshape(-1,len(index))
    result=np.array([np.bincount(index,weights=row) for row in flat])
    return result.reshape(values.shape[:-1]+(int(index.max()+1),))


def analyse(paths,out,regions_path=None,pixel_groups=None):
    out=Path(out)
    if out.exists(): raise ValueError('Use a fresh output directory')
    n,b,groups,identity=read_templates(paths)
    if pixel_groups is not None:
        if regions_path: raise ValueError('Region masks must be supplied in aggregated coordinates separately')
        n=aggregate_pixels(n,pixel_groups); b=aggregate_pixels(b,pixel_groups)
        for items in groups.values():
            for item in items: item['signal']=aggregate_pixels(item['signal'],pixel_groups)
    regions={}
    if regions_path:
        with np.load(regions_path,allow_pickle=False) as d: regions={k:d[k] for k in d.files}
    rows=[]; diagnostics=[]
    for (duration,rate,scenario),items in sorted(groups.items()):
        # States must align: common response spec/RNG/calibration states, independent dynamics seeds.
        states=combined_states(items)
        beta,nll=fixed_rate_fit(n,states,b)
        row=dict(duration_myr=duration,rate_msun_myr=rate,scenario=scenario,background_scale=beta,
            negative_log_likelihood=nll,expected_signal_stars=float(states.mean(axis=0).sum()),signal_amplitude=1.)
        rows.append(row)
        diagnostics.append(dict(duration_myr=duration,rate_msun_myr=rate,scenario=scenario,
            noise=noise_report([i['signal'].mean(axis=0) for i in items],states.mean(axis=0)+beta*b,regions),
            inputs=[{k:i[k] for k in ('path','sha256')} for i in items]))
    frame=pd.DataFrame(rows); best=frame.negative_log_likelihood.min()
    frame['q']=2*(frame.negative_log_likelihood-best)
    rates=frame.groupby('rate_msun_myr',as_index=False).q.min()
    times=frame.groupby('duration_myr',as_index=False).q.min()
    beta0,nll0=fixed_rate_fit(n,np.zeros((1,len(n))),b)
    out.mkdir(parents=True)
    frame.to_csv(out/'joint_grid.csv',index=False); rates.to_csv(out/'rate_profile.csv',index=False); times.to_csv(out/'time_profile.csv',index=False)
    report=dict(model_version=identity[-1],reporting_level=2,physical_rate_measurement_validated=False,
        pixel_groups=None if pixel_groups is None else np.asarray(pixel_groups).tolist(),
        observation_identity=identity[0],candidates=frame.to_dict(orient='records'),diagnostics=diagnostics,
        background_only=dict(scale=beta0,negative_log_likelihood=nll0),
        limitations=['Conditional grid only; no calibrated intervals or detection significance.',
            'T plateau is not automatically a lower limit.',
            'Independent-template noise diagnostic is not a substitute for refinement and injection recovery.',
            'Population/profile ranges must be declared; no automatic population-product interpretation.'])
    (out/'joint_grid.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).parent/'outputs/.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    prof=frame.groupby(['duration_myr','rate_msun_myr'],as_index=False).q.min()
    im=axes[0].scatter(prof.duration_myr/1000,prof.rate_msun_myr,c=prof.q,cmap='viridis_r')
    fig.colorbar(im,ax=axes[0],label='Relative -2 log likelihood')
    axes[0].set(xlabel='Lookback window (Gyr)',ylabel='Escape rate (solar masses/Myr)',title='Sampled grid; population profiled')
    axes[1].plot(rates.rate_msun_myr,rates.q,'o-'); axes[1].set(xlabel='Escape rate',ylabel='q',title='Profile over T and population')
    axes[2].plot(times.duration_myr/1000,times.q,'o-'); axes[2].set(xlabel='Lookback window (Gyr)',ylabel='q',title='Profile over rate and population')
    fig.suptitle('Conditional discrete grid — no calibrated confidence contours')
    fig.savefig(out/'joint_grid.png',dpi=150); plt.close(fig)
    from grid_plot import plot_grid
    plot_grid(frame.to_dict(orient='records'),out/'heatmap.png')
    (out/'REPORT.md').write_text('# Conditional T–MLR grid\n\n![Grid](joint_grid.png)\n\nReporting level 2: conditional assumptions; not a validated physical rate measurement.\n\nSee joint_grid.csv, profiles and joint_grid.json for diagnostics and limitations. No confidence intervals or T lower bound are inferred automatically.\n',encoding='utf-8')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--templates',type=Path,nargs='+',required=True); p.add_argument('--out',type=Path,required=True)
    p.add_argument('--regions',type=Path,help='Predeclared boolean masks in fitted-pixel order, NPZ')
    a=p.parse_args(); analyse(a.templates,a.out,a.regions)
