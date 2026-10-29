"""Reproducible sensitivity, numerical-resolution and conditional fit checks."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from importlib.metadata import version
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

from provenance import sha256_file

HERE=Path(__file__).resolve().parent


def command(script,args,log):
    log=Path(log); log.parent.mkdir(parents=True,exist_ok=True)
    with log.open('w',encoding='utf-8') as f:
        r=subprocess.run([sys.executable,str(HERE/script),*map(str,args)],cwd=HERE,stdout=f,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
    if r.returncode:
        raise RuntimeError(f'{script}: '+ '\n'.join(log.read_text(encoding='utf-8',errors='replace').splitlines()[-8:]))


def absolute_spec(path):
    path=Path(path); spec=json.loads(path.read_text(encoding='utf-8')); root=path.resolve().parent
    for k in ('population','quality','cache_dir','rate_history'):
        if spec.get(k): spec[k]=str((root/spec[k]).resolve())
    if spec['completeness'].get('path'): spec['completeness']['path']=str((root/spec['completeness']['path']).resolve())
    for k in ('path','metadata'): spec['errors'][k]=str((root/spec['errors'][k]).resolve())
    return spec


def check_response(path):
    path=Path(path); meta=json.loads((path/'response.json').read_text(encoding='utf-8'))
    for key,expected in meta['input_sha256'].items():
        if sha256_file(key)!=expected: raise ValueError('Response input changed: '+key)
    for key,expected in meta['source_sha256'].items():
        if sha256_file(HERE/key)!=expected: raise ValueError('Response code changed: '+key)
    for key,expected in meta.get('versions',{}).items():
        if version(key)!=expected: raise ValueError('Response dependency changed: '+key)
    survey=meta.get('completeness',{})
    if survey.get('path') and survey.get('input_sha256'):
        if sha256_file(survey['path'])!=survey['input_sha256']:
            raise ValueError('Completeness map changed')
    if sha256_file(path/'calibrated_counts.npz')!=meta['template_sha256']:
        raise ValueError('Response template checksum mismatch')
    return meta


def clear_empty_interrupted_simulation(dest,root):
    """After stopping the previous run, recover its empty output placeholder only."""
    dest,root=Path(dest),Path(root)
    if not dest.exists(): return
    if dest.resolve().parent!=root.resolve(): raise ValueError('Unexpected simulation directory')
    if any(dest.iterdir()):
        raise ValueError('Incomplete simulation contains files; preserve them and use a fresh validation directory')
    dest.rmdir()  # Nonrecursive: never delete partial scientific data.


def format_change(value):
    return 'undefined (zero baseline amplitude)' if value is None else f'{value:.2%}'


def response_job(name,spec,simulation,selection,out):
    dest=out/name; dest.mkdir(parents=True,exist_ok=True)
    response=dest/'response'
    if (response/'response.json').exists(): check_response(response)
    else: command('response_model.py',['--selection',selection,'--run',simulation,'--config',spec,'--out',response],dest/'response.log')
    if not (dest/'rate.json').exists():
        command('fit_rate.py',['--input',response/'calibrated_counts.npz','--out',dest/'rate.json'],dest/'fit.log')
    r=json.loads((dest/'rate.json').read_text(encoding='utf-8'))
    with np.load(response/'calibrated_counts.npz',allow_pickle=False) as d:
        r['reference_selected_stars']=float(d['signal'].sum())
    return dict(name=name,response=str(response.resolve()),rate=r)


def numerical_comparison(reference,alternative,tolerance):
    a,b=reference['rate'],alternative['rate']
    rate_change=abs(b['mdot_msun_myr']/a['mdot_msun_myr']-1) if a['mdot_msun_myr'] else None
    count_change=abs(b['reference_selected_stars']/a['reference_selected_stars']-1)
    with np.load(Path(reference['response'])/'calibrated_counts.npz',allow_pickle=False) as da, np.load(Path(alternative['response'])/'calibrated_counts.npz',allow_pickle=False) as db:
        sa,sb=da['signal'],db['signal']
        if not np.array_equal(da['counts'],db['counts']): raise ValueError('Convergence tests must use identical observations')
        tv=float(.5*np.sum(np.abs(sa/sa.sum()-sb/sb.sum())))
    return dict(reference=reference['name'],alternative=alternative['name'],relative_rate_change=rate_change,
        relative_selected_count_change=count_change,normalized_signal_total_variation=tv,tolerance=tolerance,
        within_rate_and_count_tolerance=bool(rate_change is not None and max(rate_change,count_change)<=tolerance),
        interpretation='A two-resolution diagnostic, not proof of asymptotic convergence or control of every spatial bin.')


def run(bundle,baseline,simulation,out,*,workers=2,particles=4000,trials=200,tolerance=.05):
    bundle,baseline,simulation,out=[Path(p).resolve() for p in (bundle,baseline,simulation,out)]
    bm=json.loads((bundle/'calibration_manifest.json').read_text(encoding='utf-8'))
    base_spec=bundle/bm['baseline']; selection=baseline/'selection'
    simmeta=json.loads((simulation/'run.json').read_text(encoding='utf-8'))
    if particles<=simmeta['config']['n_particles'] or particles%2: raise ValueError('Numerical test needs a larger even tracer count')
    if workers<1 or trials<20 or not 0<tolerance<1: raise ValueError('Invalid validation settings')
    plan=dict(bundle_manifest=sha256_file(bundle/'calibration_manifest.json'),
        scenario_configs={k:sha256_file(bundle/k) for k in bm['variants']},
        source={f:sha256_file(HERE/f) for f in ('validate_science.py','fit_diagnostics.py','response_model.py','population_model.py','selection_model.py','fit_rate.py','stream_model.py','dynamics.py')},
        original_simulation=sha256_file(simulation/'run.json'),selection=sha256_file(selection/'selection.json'),
        particles=particles,bootstrap_trials=trials,tolerance=tolerance)
    if out.exists():
        if not (out/'plan.json').exists() or json.loads((out/'plan.json').read_text(encoding='utf-8'))!=plan:
            raise ValueError('Use a fresh validation directory; existing plan differs')
    else:
        out.mkdir(parents=True); (out/'plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
    check_response(baseline/'response')
    base_rate=json.loads((baseline/'conditional_rate.json').read_text(encoding='utf-8'))
    with np.load(baseline/'response/calibrated_counts.npz') as d: base_rate['reference_selected_stars']=float(d['signal'].sum())
    reference=dict(name=Path(bm['baseline']).stem,response=str(baseline/'response'),rate=base_rate)
    models={reference['name']:reference}; configs=out/'configs'; configs.mkdir(exist_ok=True)
    mc_spec=absolute_spec(base_spec); mc_spec['observation_draws']*=2
    mc_path=configs/'double_observation_draws.json'; mc_path.write_text(json.dumps(mc_spec,indent=2),encoding='utf-8')
    cfg=dict(simmeta['config'],n_particles=particles,seed=84)
    cfg_path=configs/'more_particles.json'; cfg_path.write_text(json.dumps(cfg,indent=2),encoding='utf-8')

    def particle_job():
        dest=out/'particle_simulation'
        if not (dest/'run.json').exists():
            clear_empty_interrupted_simulation(dest,out)
            command('stream_model.py',['--config',cfg_path,'--out',dest],out/'particle_simulation.log')
        saved=json.loads((dest/'run.json').read_text(encoding='utf-8'))
        if saved['config']!=cfg: raise ValueError('Numerical simulation config mismatch')
        return response_job('more_particles',base_spec,dest,selection,out/'numerics')

    def bootstrap_job():
        dest=out/'goodness_of_fit'
        if not (dest/'diagnostics.json').exists():
            command('fit_diagnostics.py',['--input',baseline/'response/calibrated_counts.npz','--out',dest,'--trials',trials],out/'goodness_of_fit.log')
        return json.loads((dest/'diagnostics.json').read_text(encoding='utf-8'))

    numerical={}; gof=None
    with ThreadPoolExecutor(max_workers=workers) as pool:
        tasks={pool.submit(particle_job):'more_particles'}
        for name in bm['variants']:
            if name!=bm['baseline']:
                tasks[pool.submit(response_job,Path(name).stem,bundle/name,simulation,selection,out/'models')]=Path(name).stem
        tasks[pool.submit(response_job,'double_observation_draws',mc_path,simulation,selection,out/'numerics')]='double_observation_draws'
        tasks[pool.submit(bootstrap_job)]='goodness_of_fit'
        for future in as_completed(tasks):
            name=tasks[future]; result=future.result()
            if name=='goodness_of_fit': gof=result
            elif name in ('more_particles','double_observation_draws'): numerical[name]=result
            else: models[name]=result
            print('Completed validation:',name,flush=True)
    rows=[]
    for name in [Path(k).stem for k in bm['variants']]:
        m=models[name]; r=m['rate']; scenario=r['scenario']
        rows.append(dict(name=name,escape_mass_tilt=scenario['escape_mass_tilt'],parallax_zero_point_mas=scenario['parallax_zero_point_mas'],
            mdot=r['mdot_msun_myr'],lower_68=r['mdot_interval_68_asymptotic'][0],upper_68=r['mdot_interval_68_asymptotic'][1],
            nll=r['negative_log_likelihood'],reference_selected_stars=r['reference_selected_stars']))
    table=pd.DataFrame(rows); table['delta_twice_nll']=2*(table.nll-table.nll.min())
    table.to_csv(out/'sensitivity.csv',index=False)
    checks={key:numerical_comparison(reference,value,tolerance) for key,value in numerical.items()}
    report=dict(models=models,numerical_tests=checks,goodness_of_fit=gof,
        scenario_rate_range=[float(table.mdot.min()),float(table.mdot.max())],
        interpretation='Scenario ranges are assumed-model sensitivity, NOT statistical confidence intervals or a posterior for the escape mass function.',
        status='Conditional validation diagnostics only; no validated total physical rate.',
        remaining=['Population CMD/HB mismatch and evolutionary assumptions have not been resolved by this suite.',
            'Fixed background, selection and Galactic potential systematics remain.',
            'Two tracer counts with different random seeds do not establish a convergence sequence.'])
    (out/'validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    plot_sensitivity(table,out/'sensitivity.png')
    lines=['# Conditional science-model validation','',report['status'],'',
        f"Nine-scenario rate range: {table.mdot.min():.3g}–{table.mdot.max():.3g} Msun/Myr. This is NOT a confidence interval.",'',
        '![Sensitivity](sensitivity.png)','','## Numerical diagnostics','']
    for key,r in checks.items():
        lines.append(f"- {key}: rate change {format_change(r['relative_rate_change'])}; selected-count change {r['relative_selected_count_change']:.2%}; normalized-map total variation {r['normalized_signal_total_variation']:.3f}; within {tolerance:.0%} rate/count tolerance: {r['within_rate_and_count_tolerance']}.")
    lines+=['','## Conditional goodness of fit','',f"Refitted-bootstrap p = {gof['p_conditional_plus_one']:.4g}, B={trials}; background shape and selection are fixed. This is not a full-pipeline p-value.",'',
        '![Conditional bootstrap](goodness_of_fit/goodness_of_fit.png)','','## Unresolved assumptions','',*['- '+s for s in report['remaining']]]
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return report


def plot_sensitivity(table,path):
    os.environ.setdefault('MPLCONFIGDIR',str(HERE/'outputs/.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for tilt,group in table.groupby('escape_mass_tilt'):
        group=group.sort_values('parallax_zero_point_mas')
        axes[0].errorbar(group.parallax_zero_point_mas,group.mdot,
            yerr=[group.mdot-group.lower_68,group.upper_68-group.mdot],marker='o',capsize=3,label=f'Mass tilt {tilt:+g}')
        axes[1].plot(group.parallax_zero_point_mas,group.delta_twice_nll,'o-',label=f'Mass tilt {tilt:+g}')
    axes[0].set(ylabel='Conditional escape rate (solar masses / Myr)',title='Error bars: local statistical 68% intervals')
    axes[1].set(ylabel='Relative -2 log likelihood',title='Lower = better within this scenario grid')
    spread=float(table.delta_twice_nll.max())
    # An automatically zoomed 0--0.012 axis can exaggerate negligible ranking
    # differences. Preserve order-one score context without a test threshold.
    axes[1].set_ylim(-.03,max(1.05,spread*1.05))
    axes[1].text(.5,.55,f'Maximum score difference: {spread:.3g}\nLowest / highest likelihood: {np.exp(-spread/2):.3f}',
        transform=axes[1].transAxes,ha='center',fontsize=10)
    for ax in axes:
        ax.set_xlabel('Assumed scalar parallax offset (mas)'); ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle('Assumed population / zero-point sensitivity; scenario spread is NOT a confidence interval',fontsize=11)
    fig.savefig(path,dpi=160); plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('bundle','baseline','simulation','out'): p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--workers',type=int,default=2); p.add_argument('--particles',type=int,default=4000)
    p.add_argument('--gof-trials',type=int,default=200); p.add_argument('--tolerance',type=float,default=.05)
    a=p.parse_args(); r=run(a.bundle,a.baseline,a.simulation,a.out,workers=a.workers,particles=a.particles,trials=a.gof_trials,tolerance=a.tolerance)
    print(json.dumps(dict(scenario_rate_range=r['scenario_rate_range'],numerical_tests=r['numerical_tests']),indent=2))
