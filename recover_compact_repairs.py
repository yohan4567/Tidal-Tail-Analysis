"""Recover failed preparations in compact paths; never rewrite the failed version."""
import argparse,json,os,shutil,subprocess,sys,zipfile
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from provenance import sha256_file
from transfer_status_io import write_json

HERE=Path(__file__).resolve().parent
TARGETS=('NGC6341','NGC7099','NGC7078','NGC6139','NGC5904','NGC6093')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);write_json(p,d,status_file=p.name=='status.json')
def rebase(value,old,new):
    if isinstance(value,str) and value.startswith(str(old)+os.sep):return str(new)+value[len(str(old)):]
    if isinstance(value,list):return [rebase(v,old,new) for v in value]
    if isinstance(value,dict):return {k:rebase(v,old,new) for k,v in value.items()}
    return value
def semantic_spec(spec):
    spec=json.loads(json.dumps(spec))
    for k in ('population','quality','cache_dir'):spec.pop(k,None)
    for k in ('path','metadata'):spec['errors'].pop(k,None)
    return spec

def stage(previous,dest,cache,pool=False):
    previous,dest,cache=map(lambda p:Path(p).resolve(),(previous,dest,cache))
    identity=dict(source=str(previous),source_settings_sha256=sha256_file(previous/'bright_repair.json'),
                  cache_root=str(cache),spatial_pooling=pool,preparer_sha256=sha256_file(Path(__file__)))
    if dest.exists():
        old_identity=read(dest/'relocation.json')['identity'] if (dest/'relocation.json').exists() else {}
        # Controller-only bug fixes do not invalidate source-frozen numerical inputs.
        compatible=lambda value:{k:v for k,v in value.items() if k!='preparer_sha256'}
        if compatible(old_identity)!=compatible(identity):
            raise ValueError('Existing compact preparation has a different identity; use a fresh root')
        return
    settings=read(previous/'bright_repair.json');cid=settings['cluster']
    shutil.copytree(previous,dest,ignore=shutil.ignore_patterns('cache','response_cache','preflight',
        'empirical_repair_v1','__pycache__','repair_manifest.json','input_status.json','preflight.json','*.log'))
    settings=rebase(settings,previous,dest)
    bundle=Path(settings['reference']).parent;validation=None
    if pool:
        import pandas as pd
        from selection_model import QualityMap
        from quality_repair import validate_pooling
        quality,validation=validate_pooling(QualityMap.load(bundle/'quality'),pd.read_csv(bundle/'quality_holdout_cells.csv'))
        if not validation['screening_passed']:raise ValueError('Spatial-pooling held-out check failed')
        quality.save(bundle/'quality_pooled');dump(dest/'pooling_validation.json',validation)
        code=Path(settings['snapshot']);p=code/'selection_model.py';text=p.read_text(encoding='utf-8')
        needle='        return cls(arrays,meta,draw_seed)'
        if text.count(needle)!=1:raise ValueError('Unrecognized quality-map loader')
        text=text.replace(needle,"        if meta.get('quality_pooling'):\n            from quality_repair import SpatialPooledQualityMap\n            return SpatialPooledQualityMap(arrays,meta,draw_seed)\n"+needle)
        p.write_text(text,encoding='utf-8')
        for name in ('response_model.py','run_joint_study.py'):
            p=code/name;text=p.read_text(encoding='utf-8')
            if "'quality_repair.py'" not in text:text=text.replace("'selection_model.py',","'selection_model.py','quality_repair.py',",1)
            p.write_text(text,encoding='utf-8')
    for p in bundle.glob('tilt_*.json'):
        before=read(p);after=rebase(before,previous,dest)
        after['cache_dir']=str(cache/cid)
        if pool:after['quality']='quality_pooled'
        if semantic_spec(before)!=semantic_spec(after):raise ValueError('Relocation changed scientific observation settings')
        worst=Path(after['cache_dir'])/('f'*64)/'age_contributions.csv'
        if os.name=='nt' and len(str(worst))>=248:raise ValueError('Compact cache still exceeds conservative Windows path budget')
        dump(p,after)
    dump(dest/'bright_repair.json',settings)
    dump(dest/'relocation.json',dict(identity=identity,source_code_sha256={p.name:sha256_file(p) for p in (previous/'frozen_code').glob('*.py')},
        policy='Runtime paths only, except the explicit held-out-tested M80 spatial-quality pooling; donor feature distances and count thresholds unchanged.',
        pooling_validation_passed=validation['screening_passed'] if validation else None))

def finish(dest):
    dest=Path(dest).resolve();settings=read(dest/'bright_repair.json');code=Path(settings['snapshot'])
    os.environ['GC_TIDAL_STUDY']=settings['study'];sys.path.insert(0,str(code))
    from empirical_population import prepare
    from population_model import Population
    from selection_model import QualityMap
    from response_group import run_group
    from validate_science import command,check_response
    cid=settings['cluster'];source=Path(settings['source'])/cid;inp=Path(settings['calibration_inputs_root'])
    shutil.copy2(Path(settings['parent_repair'])/'mf_source.html',dest/'mf_source.html')
    reports=[]
    for offset in (0.,-.5,-1.):
        folder=dest/cid/'empirical_repair_v1'/('offset_'+f'{offset:g}'.replace('-','minus').replace('.','p'))
        prior_path=folder/'population_bundle/population_prior.json'
        prior=read(prior_path) if prior_path.exists() else prepare(cid,inp/'isochrone.dat',settings['reference'],dest/'mf_source.html',prior_path.parent,offset=offset,order=5)
        if prior['reference_config_sha256']!=sha256_file(settings['reference']):raise ValueError('Population reference changed')
        for cfg in prior['response_configs']:
            pop=Population.load(read(cfg)['population'])
            if not (pop.table.number_per_msun>0).all() or not np.isclose(np.sum(pop.table.number_per_msun*pop.table.current_mass_msun),1.,rtol=1e-12):raise ValueError('Invalid population normalization')
        reports.append(dict(offset=offset,prior=str(prior_path),nodes=5))
    c=read(source/'coupled_dynamics.json');c.update(duration_myr=2.,n_particles=16,release_sampling='stratified',python_integrator='cartesian_batch',batch_size=4,integration_rtol=1e-11,integration_atol=1e-13)
    smoke=dest/'preflight';smoke.mkdir(exist_ok=True)
    if not (smoke/'simulation/run.json').exists():
        dump(smoke/'config.json',c);command('stream_model.py',['--config',smoke/'config.json','--out',smoke/'simulation'],smoke/'integration.log')
    elif read(smoke/'simulation/run.json')['config']!=c:raise ValueError('Existing small simulation config differs')
    middle=read(reports[0]['prior'])['response_configs'][2]
    run_group(source/'selection',smoke/'simulation',[middle],[smoke/'response'],study=Path(settings['study']))
    response=check_response(smoke/'response')
    old=Path(read(dest/'relocation.json')['identity']['source'])/'preflight/response/calibrated_counts.npz'
    same=None
    if old.exists() and not (dest/'pooling_validation.json').exists():
        with np.load(old) as a,np.load(smoke/'response/calibrated_counts.npz') as b:
            same=all(np.array_equal(a[k],b[k]) for k in ('counts','background','signal','signal_draws'))
        if not same:raise ValueError('Path-only repair changed the actual numerical response')
    quality=QualityMap.load(Path(settings['reference']).parent/read(settings['reference'])['quality'])
    preflight=dict(stage='passed_small_real_input_response',tracers=16,duration_myr=2.,
        response_template_sha256=response['template_sha256'],path_only_response_identical=same,
        physical_rate_measurement_validated=False,scope='Small real-input check; full-grid and physical validation remain.')
    dump(dest/'preflight.json',preflight)
    state=dict(stage='ready_for_conditional_inference',reference=settings['reference'],population_checks=reports,
        calibration_inputs_root=str(inp),representative_response_passed=True,preflight=preflight,
        physical_rate_measurement_validated=False)
    dump(dest/'input_status.json',{cid:state})
    files=[*code.glob('*.py'),*list((code/'configs').rglob('*.json')),
        *[p for p in (dest/cid).rglob('*') if p.is_file()],dest/'mf_source.html',dest/'relocation.json',
        source/'coupled_dynamics.json',source/'selection/selection.json',source/'selection/maps.npz',Path(settings['study'])/'step3_tail_search.py']
    manifest=dict(source=settings['source'],snapshot=str(code),study=settings['study'],clusters=[cid],
        created_utc=datetime.now(timezone.utc).isoformat(),physical_rate_measurement_validated=False,
        policy='compact_paths_v2; conditional, explicitly validated spatial pooling only where declared',
        frozen_sha256={str(p.resolve()):sha256_file(p) for p in files})
    dump(dest/'repair_manifest.json',manifest);print(json.dumps(preflight),flush=True)

def run(root,previous,clusters):
    root,previous=map(lambda p:Path(p).resolve(),(root,previous));root.mkdir(parents=True,exist_ok=True)
    states=read(root/'status.json').get('clusters',{}) if (root/'status.json').exists() else {}
    def update(cid,stage,**extra):
        values=dict(stage=stage,updated_utc=datetime.now(timezone.utc).isoformat(),physical_rate_measurement_validated=False)
        values.update(extra)
        if stage!='repair_failed':states.setdefault(cid,{}).pop('error',None)
        states.setdefault(cid,{}).update(values)
        dump(root/'status.json',dict(clusters=states));print(cid,stage,extra.get('error',''),flush=True)
    def command(script,args,cid):
        from run_deadline_queue import memory_snapshot
        while (m:=memory_snapshot()) and m['available_gib']<8:
            update(cid,'waiting_for_RAM',available_gib=m['available_gib']);__import__('time').sleep(30)
        with (root/(cid+'.log')).open('a',encoding='utf-8') as log:
            r=subprocess.run([sys.executable,'-u',str(script),*map(str,args)],cwd=HERE,stdout=log,stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
        if r.returncode:raise RuntimeError('\n'.join((root/(cid+'.log')).read_text(encoding='utf-8',errors='replace').splitlines()[-10:]))
    for cid in clusters:
        try:
            dest=root/cid;old=previous/'model_repairs'/cid/'attempt_001'
            update(cid,'staging_compact_repair');stage(old,dest,HERE/'outputs/rc2',pool=cid=='NGC6093')
            if not (dest/'repair_manifest.json').exists():
                update(cid,'testing_small_response');command(Path(__file__),['--finish',dest],cid)
            settings=read(dest/'bright_repair.json');source=Path(settings['source'])/cid
            state=read(dest/'input_status.json')[cid]
            job=dict(cluster=cid,code=settings['snapshot'],study=settings['study'],config=str(source/'coupled_dynamics.json'),
                selection=str(source/'selection'),isochrone=str(Path(settings['calibration_inputs_root'])/'isochrone.dat'),reference=state['reference'],
                html=str(dest/'mf_source.html'),reuse_grid=None,cluster_age_myr=read(Path(settings['calibration_inputs_root'])/'isochrone.json')['age_years']/1e6,
                frozen_sha256=read(dest/'repair_manifest.json')['frozen_sha256'],target_parameter_posterior_computed=False,physical_rate_measurement_validated=False)
            if cid in ('NGC5904','NGC6093'):
                reuse=HERE/'outputs/cluster_survey_29_v1'/cid/'empirical_direct_v1/offset_0/grid'
                if not (reuse/'plan.json').exists():raise ValueError('Existing full-duration dynamics required')
                job['reuse_grid']=str(reuse)
            jobpath=dest/'job.json'
            if jobpath.exists() and read(jobpath)!=job:raise ValueError('Existing compact job changed')
            dump(jobpath,job)
            if cid in ('NGC5904','NGC6093'):
                update(cid,'testing_full_duration_support');command(HERE/'run_priority_model.py',['--job',jobpath,'--out',dest/'support','--workers',1,'--validate-only'],cid)
            package=root/'pkgs'/(cid+'_v2')
            if not package.with_suffix('.zip').exists():
                update(cid,'exporting_package');command(HERE/'portable_cluster.py',['export','--repair',dest,'--cluster',cid,'--out',package,'--analysis-kind','model_comparison'],cid)
            validation=root/'v'/cid
            if not (validation/'validation.json').exists():
                update(cid,'checking_relocated_response');command(HERE/'validate_transfer_package.py',['--bundle',package,'--simulation',dest/'preflight/simulation','--out',validation],cid)
            if read(validation/'validation.json')['stage']!='relocated_small_response_passed':raise ValueError('Relocated test incomplete')
            archive=package.with_suffix('.zip')
            with zipfile.ZipFile(archive) as z:
                if z.testzip() is not None:raise ValueError('ZIP integrity failed')
            entry=dict(archive=str(archive),bytes=archive.stat().st_size,archive_sha256=sha256_file(archive),zip_crc_verified=True,
                validation=str(validation/'validation.json'),analysis_kind='model_comparison',job=str(jobpath),physical_rate_measurement_validated=False)
            index=read(root/'READY_PACKAGES.json') if (root/'READY_PACKAGES.json').exists() else dict(ready_packages={})
            index['ready_packages'][cid]=entry;dump(root/'READY_PACKAGES.json',index)
            if cid not in ('NGC5904','NGC6093'):
                globalpath=HERE/'outputs/other_pc_all_remaining_20261101_v1/READY_PACKAGES.json';globalindex=read(globalpath)
                globalindex['ready_packages'][cid]=entry;globalindex['updated_utc']=datetime.now(timezone.utc).isoformat();dump(globalpath,globalindex)
            update(cid,'conditional_model_package_ready',**entry)
        except Exception as exc:update(cid,'repair_failed',error=str(exc))
    dump(root/'finished.json',dict(clusters=states,finished_utc=datetime.now(timezone.utc).isoformat()))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path);p.add_argument('--previous',type=Path)
    p.add_argument('--clusters',nargs='+',choices=TARGETS,default=list(TARGETS));p.add_argument('--finish',type=Path);p.add_argument('--run-local',action='store_true');a=p.parse_args()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    if a.finish:finish(a.finish)
    else:
        if a.root is None or a.previous is None:p.error('--root and --previous required')
        from repair_seven_clusters import worker_lock
        with worker_lock(a.root):run(a.root,a.previous,a.clusters)
        if a.run_local:
            subprocess.run([sys.executable,'-u',str(HERE/'dispatch_ready_local_repairs.py'),'--root',str(a.root),'--workers','12'],cwd=HERE,check=True)
