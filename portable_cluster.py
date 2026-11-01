"""Export verified prepared inputs, then bind and run them on another computer."""
import argparse,hashlib,json,os,re,shutil,sys,zipfile
from datetime import datetime,timezone
from importlib.metadata import version
from pathlib import Path


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def dump(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2),encoding='utf-8')


def export(repair,cluster,out):
    repair,out=Path(repair).resolve(),Path(out).resolve();info=read(repair/'repair_manifest.json')
    if out.exists() or out.with_suffix('.zip').exists():raise FileExistsError('Use a fresh export path')
    for p,sha in info['frozen_sha256'].items():
        if digest(p)!=sha:raise ValueError('Frozen source changed: '+p)
    state=read(repair/'input_status.json')[cluster]
    if state['stage']!='ready_for_conditional_inference':raise ValueError('Verified prepared inputs required')
    repo=Path(__file__).resolve().parent;source=Path(info['source'])/cluster;out.mkdir(parents=True)
    def copy(p,q):
        q=out/q;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
    code=out/'code';code.mkdir()
    for p in Path(info['snapshot']).glob('*.py'):shutil.copy2(p,code/p.name)
    shutil.copytree(Path(info['snapshot'])/'configs',code/'configs')
    rc=Path(info['snapshot'])/'.galpyrc'
    if rc.exists():shutil.copy2(rc,code/'.galpyrc')
    shutil.copytree(source/'selection',out/'inputs/selection')
    shutil.copytree(Path(state.get('calibration_inputs_root',source/'calibration_inputs')),out/'inputs/calibration_inputs',ignore=shutil.ignore_patterns('*.parts','*.log','*.download.vot'))
    shutil.copytree(Path(state['reference']).parent,out/'inputs/calibration')
    copy(source/'coupled_dynamics.json',Path('inputs/coupled_dynamics.json'))
    for item in state['population_checks']:
        prior=Path(item['prior']);shutil.copytree(prior.parent,out/'inputs/populations'/prior.parent.parent.name)
    copy(repair/'mf_source.html',Path('inputs/mf_source.html'))
    reference=read(state['reference']);completeness=Path(reference['completeness']['path'])
    copy(completeness,Path('inputs/completeness/allsky_uniq_10.fits'))
    if completeness.with_suffix(completeness.suffix+'.json').exists():
        copy(completeness.with_suffix(completeness.suffix+'.json'),Path('inputs/completeness/allsky_uniq_10.fits.json'))
    study=Path(info['study'])
    for name in ('step3_tail_search.py','skymath.py','catalogs/master_clusters.csv'):
        copy(study/name,Path('external_study')/name)
    selection=read(source/'selection/selection.json');field=Path(selection['source_field'])
    if digest(field)!=selection['source_field_sha256']:raise ValueError('Observed source field changed')
    for suffix in ('_field.csv.gz','_field.json','_missing_tiles.json'):
        p=field.parent/(cluster+suffix)
        if p.exists():copy(p,Path('external_study/gaia_fields')/p.name)
    copy(Path(__file__),Path('portable_cluster.py'));copy(repo/'LICENSE',Path('LICENSE'))
    guide=repo/'docs/PORTABLE_NGC362.md'
    if guide.exists():copy(guide,Path('RUN_GUIDE.md'))
    copy(repair/'repair_manifest.json',Path('provenance/source_repair_manifest.json'))
    if (repair/'preflight.json').exists():copy(repair/'preflight.json',Path('provenance/representative_preflight.json'))
    packages=('galpy','astropy','numpy','scipy','matplotlib','pandas','astropy-healpix','h5py')
    (out/'requirements.txt').write_text('\n'.join(f'{p}=={version(p)}' for p in packages)+'\n')
    rules=read(code/'configs/survey_populations_v1.json')['clusters']
    extra=read(code/'configs/population_supplements_v1.json')['clusters'].get(cluster,{})
    age=dict(rules.get(cluster,{}),**extra)['age_gyr']*1000.
    files={str(p.relative_to(out)).replace('\\','/'):digest(p) for p in out.rglob('*') if p.is_file()}
    manifest=dict(schema_version=1,cluster=cluster,cluster_age_myr=age,created_utc=datetime.now(timezone.utc).isoformat(),
        source_repair_manifest_sha256=digest(repair/'repair_manifest.json'),files_sha256=files,
        original_selection_sha256=selection['selection_file_sha256'],original_field_sha256=selection['source_field_sha256'],
        quality_subdirectory=Path(reference['quality']).name,
        representative_response_passed=bool(state.get('representative_response_passed',False)),
        scope='Prepared observations and conditional population/calibration inputs; no production simulation outputs or validated absolute physical rate included.')
    dump(out/'portable_manifest.json',manifest)
    with zipfile.ZipFile(out.with_suffix('.zip'),'x',compression=zipfile.ZIP_DEFLATED,compresslevel=1) as z:
        for p in sorted(out.rglob('*')):
            if p.is_file():z.write(p,Path(out.name)/p.relative_to(out))
    result=dict(directory=str(out),archive=str(out.with_suffix('.zip')),bytes=out.with_suffix('.zip').stat().st_size,
        archive_sha256=digest(out.with_suffix('.zip')))
    dump(out.parent/(out.name+'_export.json'),result);print(json.dumps(result,indent=2))


def verify(base):
    info=read(base/'portable_manifest.json')
    for name,sha in info['files_sha256'].items():
        path=(base/name).resolve()
        if not path.is_relative_to(base.resolve()) or digest(path)!=sha:raise ValueError('Package file missing or changed: '+name)
    if digest(base/'external_study/step3_tail_search.py')!=info['original_selection_sha256']:
        raise ValueError('Selection implementation changed')
    os.environ['GC_TIDAL_STUDY']=str((base/'external_study').resolve())
    sys.path.insert(0,str(base/'code'))
    from population_model import Population
    import numpy as np
    support=None;count=0
    for prior in sorted((base/'inputs/populations').glob('offset_*/population_prior.json')):
        for folder in sorted(prior.parent.glob('population_*')):
            if not folder.is_dir():continue
            p=Population.load(folder);t=p.table
            columns=[c for c in t.columns if c!='number_per_msun']
            if support is None:support=t[columns]
            elif not support.equals(t[columns]):raise ValueError('Population photometric supports differ')
            if not (t.number_per_msun>0).all() or not np.isclose(np.sum(t.number_per_msun*t.current_mass_msun),1.,rtol=1e-12):
                raise ValueError('Population support/mass normalization failed')
            count+=1
    if count!=15:raise ValueError('Expected three offsets with five population nodes each')
    print(f"Verified {len(info['files_sha256'])} files and {count} population nodes for {info['cluster']}",flush=True)
    return info


def prepare(base,out,info):
    out=out.resolve();identity=dict(package_manifest_sha256=digest(base/'portable_manifest.json'),package_root=str(base.resolve()),
        output_root=str(out),runner_sha256=digest(base/'portable_cluster.py'),cluster=info['cluster'])
    binding=out/'portable_binding.json'
    if binding.exists():
        if read(binding)!=identity:raise ValueError('Moved or changed bound run; use a new output directory')
        return
    if out.exists():raise ValueError('Choose a fresh output directory')
    out.mkdir(parents=True)
    for original in sorted((base/'inputs/populations').glob('offset_*/population_prior.json')):
        prior=read(original);dest=out/original.parent.name;configs=[]
        for old in sorted(original.parent.glob('response_*.json')):
            s=read(old);node='population_'+old.stem.split('_')[-1]
            s.update(population=str((old.parent/node).resolve()),quality=str((base/'inputs/calibration'/info.get('quality_subdirectory','quality')).resolve()),
                cache_dir=str(out/'response_cache'/original.parent.name))
            s['errors'].update(path=str((base/'inputs/calibration/parent_spatial_calibration.csv.gz').resolve()),
                metadata=str((base/'inputs/calibration/parent_spatial_calibration.csv.json').resolve()))
            s['completeness']['path']=str((base/'inputs/completeness/allsky_uniq_10.fits').resolve())
            path=dest/'configs'/old.name;dump(path,s);configs.append(str(path))
        prior.update(response_configs=configs,portable_source_prior_sha256=digest(original),
            portable_binding=identity,portable_relocation='Only runtime file locations rebound; original weighted populations and physical observation settings unchanged.')
        dump(dest/'population_prior.json',prior)
    dump(binding,identity)
    print('Prepared runtime paths:',out,flush=True)


def run(base,out,info,workers,response_workers=None):
    from run_joint_study import run as grid
    from bayesian_grid import run as infer
    from production_study import RATES,TIMES,SEEDS
    response_workers=workers if response_workers is None else response_workers
    if workers<1 or response_workers<1:raise ValueError('Positive worker counts required')
    reuse=None;status=out/'run_status.json'
    execution=dict(started_utc=datetime.now(timezone.utc).isoformat(),workers=workers,response_workers=response_workers,
        launcher=str(Path(__file__).resolve()),launcher_sha256=digest(Path(__file__)),
        package_root=str(base),package_manifest_sha256=digest(base/'portable_manifest.json'))
    dump(out/'execution_settings.json',execution)
    try:
        for offset in ('offset_0','offset_minus0p5','offset_minus1'):
            dest=out/offset;prior_path=dest/'population_prior.json';prior=read(prior_path)
            dump(status,dict(stage='running_conditional_grid',cluster=info['cluster'],offset=offset,workers=workers,response_workers=response_workers,
                updated_utc=datetime.now(timezone.utc).isoformat(),physical_rate_measurement_validated=False))
            grid(base/'inputs/coupled_dynamics.json',prior['response_configs'],base/'inputs/selection',dest/'grid',rates=RATES,times=TIMES,
                seeds=SEEDS,interval_myr=.25,workers=workers,response_workers=response_workers,cluster_age_myr=info['cluster_age_myr'],
                pixel_factor=4,integrator='cartesian_batch',batch_size=128,reuse_grid=reuse)
            infer(dest/'grid',dest/'posterior',prior_path)
            if offset=='offset_0':reuse=dest/'grid'
        dump(status,dict(stage='completed_conditional_posteriors_validation_pending',cluster=info['cluster'],
            updated_utc=datetime.now(timezone.utc).isoformat(),physical_rate_measurement_validated=False))
    except Exception as exc:
        dump(status,dict(stage='failed',error=str(exc),cluster=info['cluster'],updated_utc=datetime.now(timezone.utc).isoformat(),physical_rate_measurement_validated=False));raise


def status(out):
    try:state=read(out/'run_status.json')
    except FileNotFoundError:state={'stage':'not_started'}
    except json.JSONDecodeError:state={'stage':'status_update_in_progress'}
    result=dict(state=state,grids=[],note='Stage percentages are work counts, not overall elapsed-time percentages. Particle progress uses the latest log output, which can be buffered.')
    for plan in sorted(out.glob('offset_*/grid/plan.json')):
        p=read(plan);grid=plan.parent;total=len(p['rates'])*len(p['seeds']);maximum=max(p['times_myr'])
        particles=2*round(maximum/p['release_interval_myr']);integrated=complete=0
        for rate in p['rates']:
            for seed in p['seeds']:
                model=grid/f'rate_{rate:g}_seed_{seed}';done=(model/'simulation/run.json').exists()
                count=particles if done else 0;complete+=int(done);log=model/'simulation.log'
                if not done and log.exists():
                    with log.open('rb') as f:
                        f.seek(max(0,log.stat().st_size-8192));tail=f.read().decode('utf-8',errors='replace')
                    matches=re.findall(r'(?:Integrated |Resuming integrated tracers: )(\d+)/(\d+)',tail)
                    if matches and int(matches[-1][1])==particles:count=min(particles,int(matches[-1][0]))
                integrated+=count
        full=len(list(grid.glob(f'rate_*/population_*_T_{maximum:g}/response.json')))
        windows=len(list(grid.glob('rate_*/population_*/response.json')))
        expected=total*len(p['specs']);all_windows=expected*len(p['times_myr'])
        nodes=grid.parent/'posterior/integrated_nodes.jsonl';node_count=0
        if nodes.exists():
            with nodes.open('rb') as f:node_count=sum(line.endswith(b'\n') for line in f)
        node_total=len(p['rates'])*len(p['times_myr'])*len(p['specs'])
        percent=lambda n,d:round(100*n/d,2) if d else None
        result['grids'].append(dict(offset=grid.parent.name,dynamics_complete=complete,dynamics_total=total,
            dynamics_percent=percent(complete,total),integrated_particles=integrated,particles_total=particles*total,
            particle_percent=percent(integrated,particles*total),full_templates=full,full_templates_total=expected,
            full_templates_percent=percent(full,expected),window_templates=windows,window_templates_total=all_windows,
            window_templates_percent=percent(windows,all_windows),background_integrals_complete=node_count,
            background_integrals_total=node_total,background_integrals_percent=percent(node_count,node_total),
            posterior_ready=(grid.parent/'posterior/posterior_summary.json').exists()))
    print(json.dumps(result,indent=2))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('export','verify','prepare','run','status'))
    p.add_argument('--repair',type=Path);p.add_argument('--cluster');p.add_argument('--out',type=Path)
    p.add_argument('--workers',type=int,default=8)
    p.add_argument('--response-workers',type=int,help='Concurrent observation jobs; defaults to min(workers,8)')
    p.add_argument('--bundle',type=Path,help='Prepared package directory; defaults to this script directory')
    a=p.parse_args();base=(a.bundle or Path(__file__).resolve().parent).resolve()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    if a.workers<1 or (a.response_workers is not None and a.response_workers<1):p.error('Worker counts must be positive')
    if a.action=='export':
        if a.repair is None or a.cluster is None or a.out is None:p.error('export requires --repair, --cluster and --out')
        export(a.repair,a.cluster,a.out);return
    out=(a.out or base/'results').resolve()
    if a.action=='status':status(out);return
    info=verify(base)
    if a.action=='verify':return
    if a.action=='run' and not info.get('representative_response_passed'):
        raise ValueError('This export has not passed a representative observation response; repair its calibration before starting the production grid')
    prepare(base,out,info)
    if a.action=='run':run(base,out,info,a.workers,a.response_workers if a.response_workers is not None else min(a.workers,8))


if __name__=='__main__':main()
