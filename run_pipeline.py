"""One command for observed tail profiles, optional simulation and calibrated rate fitting."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import platform
import time
from importlib.metadata import version

from selection_bridge import DEFAULT_STUDY, load_study
from provenance import sha256_file

HERE=Path(__file__).resolve().parent


def execute(script,arguments,out,stages):
    log=out/(script.replace('.py','')+'.log')
    print(f'Running {script} ...',flush=True)
    started=time.perf_counter()
    with log.open('w',encoding='utf-8') as f:
        r=subprocess.run([sys.executable,str(HERE/script),*map(str,arguments)],cwd=HERE,stdout=f,stderr=subprocess.STDOUT)
    stages.append(dict(script=script,returncode=r.returncode,log=log.name,elapsed_seconds=time.perf_counter()-started))
    if r.returncode:
        tail='\n'.join(log.read_text(encoding='utf-8',errors='replace').splitlines()[-8:])
        raise RuntimeError(f'{script} failed. See {log}\n{tail}')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cluster',default='Pal5')
    p.add_argument('--study',type=Path,default=DEFAULT_STUDY)
    p.add_argument('--legacy',action='store_true')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--mode',choices=('observed','calibrated','weighted'),default='observed')
    p.add_argument('--response-config',type=Path,help='Hybrid v2 population/completeness/quality/error configuration')
    sim=p.add_mutually_exclusive_group()
    sim.add_argument('--simulate',action='store_true')
    sim.add_argument('--reuse-run',type=Path)
    p.add_argument('--duration-myr',type=float,default=1000)
    p.add_argument('--particles',type=int,default=2000)
    p.add_argument('--mdot',type=float,default=1.)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--population',type=Path)
    p.add_argument('--completeness',type=Path)
    p.add_argument('--a-g',type=float)
    p.add_argument('--reddening',type=float)
    p.add_argument('--calibration-description')
    p.add_argument('--synthetic-calibration',action='store_true')
    p.add_argument('--null-trials',type=int,default=0)
    p.add_argument('--bin-deg',type=float)
    p.add_argument('--half-width-deg',type=float)
    p.add_argument('--alpha',type=float,default=.05)
    p.add_argument('--min-coverage',type=float,default=.7)
    p.add_argument('--min-consecutive',type=int,default=2)
    a=p.parse_args()
    a.study=a.study.resolve(); a.out=a.out.resolve()
    for key in ('population','completeness','reuse_run','response_config'):
        value=getattr(a,key)
        if value is not None:
            setattr(a,key,value.resolve())
    if a.mode=='weighted':
        if not (a.simulate or a.reuse_run) or a.response_config is None:
            p.error('Weighted mode needs --simulate/--reuse-run and --response-config')
        if any(getattr(a,k) is not None for k in ('population','completeness','a_g','reddening','calibration_description')):
            p.error('Put population, completeness and extinction in the weighted response config')
    elif a.response_config is not None:
        p.error('--response-config requires --mode weighted')
    elif a.mode=='calibrated':
        if not (a.simulate or a.reuse_run) or any(getattr(a,k) is None for k in ('population','completeness','a_g','reddening','calibration_description')):
            p.error('Calibrated mode needs --simulate/--reuse-run, --population, --completeness, --a-g, --reddening and --calibration-description')
    elif a.synthetic_calibration or any(getattr(a,k) is not None for k in ('population','completeness','a_g','reddening','calibration_description')):
        p.error('Calibration arguments require --mode calibrated')
    T=load_study(a.study,a.legacy)
    field=T.field_file(a.cluster,'_field.csv.gz')
    if not field.exists():
        p.error(f'Planned field not available: {field}. Use --legacy only if intentionally analysing the old field.')
    a.out.mkdir(parents=True,exist_ok=False)
    stages=[]
    manifest=dict(cluster=a.cluster,mode=a.mode,status='running',stages=stages,
                  versions=dict(python=platform.python_version(),**{name:version(name) for name in ('numpy','scipy','pandas','matplotlib','astropy','galpy')}),
                  source_sha256={f.name:sha256_file(f) for f in HERE.glob('*.py')},
                  arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()})
    try:
        simulation=a.reuse_run
        if a.simulate:
            config=a.out/'config.json'
            execute('make_config.py',[a.cluster,'--catalog',a.study/'catalogs/master_clusters.csv','--out',config],a.out,stages)
            c=json.loads(config.read_text(encoding='utf-8'))
            c.update(duration_myr=a.duration_myr,n_particles=a.particles,reference_mdot_msun_myr=a.mdot,seed=a.seed)
            config.write_text(json.dumps(c,indent=2),encoding='utf-8')
            simulation=a.out/'simulation'
            execute('stream_model.py',['--config',config,'--out',simulation],a.out,stages)
        if simulation:
            meta=json.loads((simulation/'run.json').read_text(encoding='utf-8'))
            if meta['config']['name']!=a.cluster:
                raise ValueError('Simulation and observed cluster IDs differ')
            execute('model_extent.py',['--run',simulation,'--out',a.out/'model_extent'],a.out,stages)
        selection_args=['--cluster',a.cluster,'--study',a.study,'--out',a.out/'selection']
        if a.legacy:
            selection_args+=['--legacy']
        if a.mode=='calibrated':
            mock=a.out/'mock.csv'
            observe_args=['--run',simulation,'--population',a.population,'--completeness',a.completeness,
                    '--field',field,'--a-g',a.a_g,'--reddening',a.reddening,'--seed',a.seed,
                    '--calibration-description',a.calibration_description,'--out',mock]
            if a.synthetic_calibration:
                observe_args+=['--synthetic-calibration']
            execute('observe_particles.py',observe_args,a.out,stages)
            selection_args+=['--mock',mock,'--reference-mdot',meta['config']['reference_mdot_msun_myr'],
                             '--calibration-description',a.calibration_description]
        execute('selection_bridge.py',selection_args,a.out,stages)
        template_path=a.out/'selection/calibrated_counts.npz'
        if a.mode=='weighted':
            execute('response_model.py',['--selection',a.out/'selection','--run',simulation,
                    '--config',a.response_config,'--study',a.study,'--out',a.out/'response'],a.out,stages)
            response=json.loads((a.out/'response/response.json').read_text(encoding='utf-8'))
            a.synthetic_calibration=bool(response['synthetic_calibration'])
            manifest['synthetic_calibration']=a.synthetic_calibration
            template_path=a.out/'response/calibrated_counts.npz'
        extent_args=['--selection',a.out/'selection','--out',a.out/'tail_extent','--alpha',a.alpha,
                     '--min-coverage',a.min_coverage,'--min-consecutive',a.min_consecutive]
        if a.bin_deg is not None:
            extent_args+=['--bin-deg',a.bin_deg]
        if a.half_width_deg is not None:
            extent_args+=['--half-width-deg',a.half_width_deg]
        execute('tail_extent.py',extent_args,a.out,stages)
        if a.mode in ('calibrated','weighted'):
            execute('fit_rate.py',['--input',template_path,'--out',a.out/'conditional_rate.json',
                                   '--null-trials',a.null_trials,'--seed',a.seed],a.out,stages)
            execute('plot_rate_profile.py',['--input',template_path,'--out',a.out/'rate_profile','--label',a.cluster],a.out,stages)
        sr=json.loads((a.out/'selection/selection.json').read_text(encoding='utf-8'))
        er=json.loads((a.out/'tail_extent/extent.json').read_text(encoding='utf-8'))
        lines=[f'# {a.cluster}: completed analysis run','',f'Mode: {a.mode}. Legacy field: {a.legacy}.',
               f"Observed candidates: {sr['observation_cutflow']['cmd_and_missing_tiles']}; exact match to original selection: {sr['exact_selection_match']}.",
               '', '## Conditional projected extent', '',
               'The background is estimated from this field. These are conditional orbit-aligned candidates, not a physical tail boundary.', '']
        if a.synthetic_calibration:
            lines[2:2]=['**SYNTHETIC CALIBRATION — software integration test only. Do not interpret this run as a measured physical rate.**','']
        elif a.mode=='weighted' and response.get('scientific_status','').startswith('conditional_real_input'):
            lines[2:2]=['**REAL INPUTS, CONDITIONAL MODEL — population/CMD, calibration transfer and sampling are not yet validated for a physical mass-loss measurement.**',
                        '',f"Rate scope: {response['rate_scope']}",'']
        for direction,arm in er['directions'].items():
            extent=arm['connected_extent_deg']
            label=f'{extent:.3f} tangent deg' if extent is not None else 'no connected candidate at the chosen settings'
            lines.append(f'- {direction}: {label}; candidate segments: {len(arm["candidate_segments"])}.')
        lines+=['','![Tail profile](tail_extent/extent.png)','','![Selection maps](selection/selection.png)',
                '', 'See tail_extent/extent.json for coverage, width, search and background limitations.']
        if simulation:
            mr=json.loads((a.out/'model_extent/model_extent.json').read_text(encoding='utf-8'))
            lines+=['','## Model radial mass quantiles','',
                    'Model outputs depend on release duration and sampling. They are not observed lengths or a physical tail boundary.','']
            for name,arm in mr['directions'].items():
                if 'radial_quantile_deg' in arm:
                    lines.append(f'- {name}: 95% radial mass quantile {arm["radial_quantile_deg"]:.4g} deg; {arm["n_tracers"]} tracers; {arm["sampling_status"]}.')
            lines+=['','See model_extent/model_extent.json and the simulation run.json for assumptions and numerical checks.']
        if a.mode in ('calibrated','weighted'):
            rr=json.loads((a.out/'conditional_rate.json').read_text(encoding='utf-8'))
            lines+=['','## Conditional rate', '',f"Rate: {rr['mdot_msun_myr']:.4g} Msun/Myr; 68% asymptotic interval: {rr['mdot_interval_68_asymptotic']}.",
                    'Population and completeness were supplied externally. Statistical intervals exclude unrepresented calibration and model-shape uncertainties.']
            lines+=['','![Assumed rate versus relative Poisson fit](rate_profile/profile.png)']
            if a.mode=='weighted':
                lines+=['',f"Mass basis: {response['mass_basis']}. Shared availability/quality calibration draws: {response['calibration_draws']}.",
                        'The likelihood marginalizes the supplied shared calibration draws. Other unrepresented calibration/model systematics remain conditional.',
                        '', '![Weighted response](response/response.png)']
        else:
            lines+=['','No observed mass-loss rate was inferred. Use calibrated mode with documented population/completeness inputs.']
        (a.out/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
        manifest['status']='complete'
    except Exception as exc:
        manifest.update(status='failed',error=str(exc))
        raise
    finally:
        (a.out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(f'Completed: {a.out / "REPORT.md"}',flush=True)


if __name__=='__main__':
    main()
