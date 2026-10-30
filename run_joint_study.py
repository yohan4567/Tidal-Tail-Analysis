"""Run one Tmax dynamics per rate/seed, then window responses and a conditional grid."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import numpy as np

from joint_grid import analyse
from provenance import sha256_file
from validate_science import command,check_response

HERE=Path(__file__).resolve().parent


def run(config,specs,selection,out,*,rates,times,seeds=(42,84),interval_myr=1.,workers=2,cluster_age_myr=None,pixel_factor=1,integrator="cartesian",batch_size=32,response_workers=None):
    config,selection,out=map(lambda p:Path(p).resolve(),(config,selection,out))
    specs=[Path(p).resolve() for p in specs]
    c=json.loads(config.read_text(encoding='utf-8'))
    for name,values in [('rates',rates),('times',times)]:
        if not values or not np.isfinite(values).all() or min(values)<=0 or len(set(values))!=len(values):
            raise ValueError(f'Positive distinct {name} required')
    if not seeds or len(set(seeds))!=len(seeds) or any(int(s)!=s or s<0 for s in seeds):
        raise ValueError('Distinct nonnegative integer seeds required')
    if not np.isfinite(interval_myr) or interval_myr<=0 or workers<1: raise ValueError('Positive interval/workers required')
    if cluster_age_myr is None or not np.isfinite(cluster_age_myr) or cluster_age_myr<=0:
        raise ValueError('Explicit positive cluster age upper limit required; age uncertainty/history must be assessed separately')
    if max(times)>cluster_age_myr: raise ValueError('Lookback window exceeds supplied cluster age')
    if any(not np.isclose(t/interval_myr,round(t/interval_myr),rtol=0,atol=1e-7) for t in times):
        raise ValueError('All times must align with release interval')
    if not c.get('mass_evolution') or c['mass_evolution'].get('mode')!='linear_escape':
        raise ValueError('Explicit linear_escape mass conversion/basis/source required in dynamics config')
    if integrator not in ('cartesian','cartesian_batch'):raise ValueError('Unsupported Cartesian backend')
    from stream_model import validate
    sources={k:sha256_file(HERE/k) for k in ('run_joint_study.py','joint_grid.py','joint_fit.py','stream_model.py','dynamics.py','mass_evolution.py','response_model.py','population_model.py','selection_model.py','validate_window_reuse.py','response_windows.py')}
    plan=dict(model_version='prescribed_mass_history_v1',method_document_version=3,rates=list(rates),times_myr=list(times),
        seeds=list(seeds),release_interval_myr=interval_myr,cluster_age_myr=cluster_age_myr,
        config_sha256=sha256_file(config),specs={str(p):sha256_file(p) for p in specs},
        selection_sha256=sha256_file(selection/'selection.json'),maps_sha256=sha256_file(selection/'maps.npz'),source_sha256=sources,
        conditional_only=True,subset_validation='Pending independent short-window comparison; outputs cannot certify reuse automatically')
    if integrator=='cartesian_batch':
        sources.update({name:sha256_file(HERE/name) for name in ('batch_dynamics.py','mw2014_fast.py')})
        plan.update(integrator=integrator,batch_size=batch_size,integration_rtol=1e-11,integration_atol=1e-13)
    if not isinstance(pixel_factor,int) or pixel_factor<1: raise ValueError('Positive integer pixel factor required')
    if pixel_factor!=1: plan['pixel_factor']=pixel_factor
    if out.exists():
        if not (out/'plan.json').exists() or json.loads((out/'plan.json').read_text())!=plan:
            raise ValueError('Existing plan differs; choose a fresh output directory')
    else:
        out.mkdir(parents=True); (out/'plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
    configs=out/'configs'; configs.mkdir(exist_ok=True)
    jobs=[]; feasibility=[]
    for rate in rates:
        for seed in seeds:
            cc=dict(c,duration_myr=max(times),reference_mdot_msun_myr=rate,seed=seed,
                n_particles=2*round(max(times)/interval_myr),release_sampling='stratified',python_integrator=integrator)
            if integrator=='cartesian_batch':cc.update(batch_size=batch_size,integration_rtol=1e-11,integration_atol=1e-13)
            validate(cc)
            ratio=1+cc['mass_evolution']['dynamical_to_population_mass']*rate*max(times)/cc['mass_msun']
            feasibility.append(dict(rate=rate,seed=seed,past_to_present_mass_ratio=ratio,status='model_applicability_requires_review'))
            label=f'rate_{rate:g}_seed_{seed}'; cfg=configs/(label+'.json'); cfg.write_text(json.dumps(cc,indent=2))
            jobs.append((label,cfg,out/label/'simulation',cc))
    (out/'applicability.json').write_text(json.dumps(feasibility,indent=2))
    def dynamics_job(job):
        label,cfg,dest,cc=job; dest.parent.mkdir(exist_ok=True)
        if not (dest/'run.json').exists():
            args=['--config',cfg,'--out',dest]
            if dest.exists(): args+=['--resume']
            command('stream_model.py',args,dest.parent/'simulation.log')
        meta=json.loads((dest/'run.json').read_text())
        if meta['config']!=cc or sha256_file(dest/'particles.csv')!=meta['particles_sha256']:
            raise ValueError('Existing dynamics config/checksum mismatch')
        print('Dynamics complete:',label,flush=True)
        return label,dest
    with ThreadPoolExecutor(max_workers=workers) as pool: simulations=list(pool.map(dynamics_job,jobs))
    if min(times)<max(times) and not (out/'window_audit.json').exists():
        from validate_window_reuse import validate_window
        audit=validate_window(simulations[0][1],min(times))
        (out/'window_audit.json').write_text(json.dumps(audit,indent=2))
    response_jobs=[]
    for label,simulation in simulations:
        for i,spec in enumerate(specs):
            # Preserve relative calibration paths after saving the derived config elsewhere.
            specdata=json.loads(spec.read_text(encoding='utf-8'))
            for key in ('population','quality','cache_dir'):
                if specdata.get(key): specdata[key]=str((spec.parent/specdata[key]).resolve())
            for key in ('path','metadata'):
                specdata['errors'][key]=str((spec.parent/specdata['errors'][key]).resolve())
            completeness=specdata['completeness']
            for key in ('path','metadata'):
                if isinstance(completeness.get(key),str): completeness[key]=str((spec.parent/completeness[key]).resolve())
            if specdata.get('rate_history'): raise ValueError('External rate history not supported by constant coupled grid')
            for duration in [max(times)]:
                name=f'{label}_pop_{i}_T_{duration:g}'; cfg=configs/(name+'.json')
                cfg.write_text(json.dumps(dict(specdata,lookback_window_myr=duration),indent=2))
                response_jobs.append((simulation,cfg,out/label/f'population_{i}_T_{duration:g}'))
    def response_job(job):
        simulation,cfg,dest=job
        if not (dest/'response.json').exists():
            command('response_model.py',['--selection',selection,'--run',simulation,'--config',cfg,'--out',dest],dest.parent/(dest.name+'.log'))
        check_response(dest)
        meta=json.loads((dest/'response.json').read_text())
        if meta['input_sha256'].get(str(cfg.resolve()))!=sha256_file(cfg): raise ValueError('Stale response config')
        from response_windows import derive
        templates=[]
        for duration in times:
            window=dest if duration==max(times) else dest.parent/(dest.name.rsplit('_T_',1)[0]+f'_T_{duration:g}')
            if window!=dest and not (window/'response.json').exists(): derive(dest,window,duration)
            check_response(window)
            templates.append(window/'calibrated_counts.npz')
        return templates
    with ThreadPoolExecutor(max_workers=response_workers or workers) as pool:
        templates=[p for group in pool.map(response_job,response_jobs) for p in group]
    manifest=dict(templates=[str(p.resolve()) for p in templates])
    (out/'templates.json').write_text(json.dumps(manifest,indent=2))
    if not (out/'scores').exists():
        labels=None
        if pixel_factor!=1:
            with np.load(selection/'maps.npz',allow_pickle=False) as data:
                mask=data['fitmask']; x,y=np.nonzero(mask)
                labels=(x//pixel_factor)*int(np.ceil(mask.shape[1]/pixel_factor))+y//pixel_factor
        analyse(templates,out/'scores',pixel_groups=labels)
    print('Conditional joint study complete; validation remains required:',out,flush=True)
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('config','selection','out'): p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--response-configs',type=Path,nargs='+',required=True)
    p.add_argument('--rates',type=float,nargs='+',required=True)
    p.add_argument('--times-myr',type=float,nargs='+',required=True)
    p.add_argument('--seeds',type=int,nargs='+',default=[42,84])
    p.add_argument('--interval-myr',type=float,default=1.)
    p.add_argument('--cluster-age-myr',type=float,required=True)
    p.add_argument('--pixel-factor',type=int,default=1)
    p.add_argument('--workers',type=int,default=2)
    a=p.parse_args(); run(a.config,a.response_configs,a.selection,a.out,rates=a.rates,times=a.times_myr,
        seeds=a.seeds,interval_myr=a.interval_myr,workers=a.workers,cluster_age_myr=a.cluster_age_myr,pixel_factor=a.pixel_factor)
