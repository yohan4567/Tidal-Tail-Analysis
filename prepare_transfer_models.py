"""Build conditional model packages from a versioned observed-transfer plan."""
import argparse,json,os,shutil,subprocess,sys,time
from pathlib import Path
from datetime import datetime,timezone
from provenance import sha256_file
from transfer_status_io import write_json

HERE=Path(__file__).resolve().parent


def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,d):p=Path(p);write_json(p,d,status_file=p.name=='model_preparation_status.json')


def run(root,retry_failed=False,ready_index=None):
    root=Path(root).resolve();study=root/'study';source=HERE/'outputs/cluster_survey_29_v1';old=HERE/'outputs/input_repairs_20261101_v1'
    status_path=root/'model_preparation_status.json'
    states=read(status_path).get('clusters',{}) if status_path.exists() else {}
    while not (root/'transfer_plan.json').exists():time.sleep(10)
    clusters=read(root/'transfer_plan.json')['clusters']
    def update(cid,stage,**extra):
        states.setdefault(cid,{}).update(stage=stage,updated_utc=datetime.now(timezone.utc).isoformat(),**extra)
        dump(root/'model_preparation_status.json',dict(clusters=states));print(cid,stage,extra.get('error',''),flush=True)
    def available():
        try:return read(root/'status.json')['clusters']
        except (FileNotFoundError,json.JSONDecodeError):return {}
    def headroom(cid):
        from run_deadline_queue import memory_snapshot
        while True:
            m=memory_snapshot()
            if not m or m['available_gib']>=8:return
            update(cid,'waiting_for_RAM',available_gib=m['available_gib']);time.sleep(30)
    def command(script,args,log,cid):
        headroom(cid)
        with log.open('a',encoding='utf-8') as f:
            r=subprocess.run([sys.executable,'-u',str(script),*map(str,args)],cwd=HERE,stdout=f,stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
        if r.returncode:raise RuntimeError('\n'.join(log.read_text(encoding='utf-8',errors='replace').splitlines()[-8:]))
    ready=read(ready_index).get('ready_packages',{}) if ready_index is not None else {}
    done=set()
    for cid in clusters:
        entry=ready.get(cid,{})
        if entry.get('archive') and Path(entry['archive']).is_file() and Path(entry.get('validation','')).is_file():
            if read(entry['validation']).get('stage')=='relocated_small_response_passed':
                states.setdefault(cid,{}).pop('error',None)
                update(cid,'conditional_model_package_ready',archive=entry['archive'],completed_in_separate_version=True,
                    conditional_model_run_ready=True,physical_rate_measurement_validated=False);done.add(cid);continue
        if states.get(cid,{}).get('stage')=='model_preparation_failed' and not retry_failed:continue
        update(cid,'queued',physical_rate_measurement_validated=False)
    base=root/'model_base';code=base/'code'
    if not code.exists():
        original=read(old/'repair_manifest.json')
        for file,sha in original['frozen_sha256'].items():
            if sha256_file(file)!=sha:raise ValueError('Frozen source changed')
        shutil.copytree(original['snapshot'],code,ignore=shutil.ignore_patterns('__pycache__','outputs'))
        for name in ('bright_calibration.py','likelihood_comparison.py','quality_repair.py'):shutil.copy2(HERE/name,code/name)
        p=code/'prepare_science_calibration.py';text=p.read_text(encoding='utf-8').replace('allsky_uniq_10.fits','allsky_M10_hpx7.hdf5').replace("mode='multi'","mode='hpx7'").replace('Official GaiaUnlimited TCG multi-resolution map','Official GaiaUnlimited TCG fixed-resolution HEALPix order-7 map');p.write_text(text,encoding='utf-8')
        p=code/'selection_model.py';text=p.read_text(encoding='utf-8');needle='        return cls(arrays,meta,draw_seed)'
        text=text.replace(needle,"        if meta.get('quality_pooling'):\n            from quality_repair import SpatialPooledQualityMap\n            return SpatialPooledQualityMap(arrays,meta,draw_seed)\n"+needle);p.write_text(text,encoding='utf-8')
        for name in ('response_model.py','run_joint_study.py'):
            p=code/name;p.write_text(p.read_text(encoding='utf-8').replace("'selection_model.py',","'selection_model.py','quality_repair.py',",1),encoding='utf-8')
        # Unsupported cells remain unsupported. Validate actual model queries
        # rather than demanding support in colour bins never used by a model.
        p=code/'bright_calibration.py';text=p.read_text(encoding='utf-8')
        start=text.index("    if min(quality.arrays['n_parent'].min(),quality.arrays['n_available'].min())<quality.minimum:")
        end=text.index('    # The original MF table',start)
        text=text[:start]+text[end:]
        p.write_text(text,encoding='utf-8')
        shutil.copy2(old/'mf_source.html',base/'mf_source.html')
        files=[*code.glob('*.py'),*list((code/'configs').rglob('*.json')),base/'mf_source.html']
        dump(base/'repair_manifest.json',dict(source=str(root/'observed'),snapshot=str(code),study=str(study),
            preparation_source_sha256=sha256_file(Path(__file__)),
            frozen_sha256={str(p.resolve()):sha256_file(p) for p in files}))
    os.environ['GC_TIDAL_STUDY']=str(study)
    for cid in clusters:
        if cid in done or (states.get(cid,{}).get('stage')=='model_preparation_failed' and not retry_failed):continue
        try:
            while not available().get(cid,{}).get('observation_ready'):
                s=available().get(cid,{}).get('stage','')
                if s in ('field_missing_or_observation_failed','observation_preparation_failed'):raise ValueError('Observed field/selection failed; inspect observation logs')
                update(cid,'waiting_for_observed_field');time.sleep(30)
            dest=root/'observed'/cid;inp=dest/'calibration_inputs';inp.mkdir(exist_ok=True)
            from fetch_isochrone import validate_returned_model
            from population_model import read_isochrone
            for suffix in ('dat','json'):
                target=inp/('isochrone.'+suffix)
                if not target.exists():shutil.copy2(source/cid/'calibration_inputs'/target.name,target)
            meta=read(inp/'isochrone.json');validate_returned_model(read_isochrone(inp/'isochrone.dat'),meta['age_years'],meta['metallicity_mh'])
            cfg=read(dest/'adopted_dynamics.json');cfg['mass_evolution']=dict(mode='linear_escape',dynamical_to_population_mass=1.,population_mass_basis='isochrone_current_stellar_mass',source_description=f'{cid}: conditional living-star baseline f=1, fixed Plummer scale; no remnant or wind correction.')
            dump(dest/'coupled_dynamics.json',cfg);report=read(dest/'selection/selection.json')
            parent=inp/'parent_random10.csv.gz';parent_meta=parent.with_suffix('.json')
            existing=old/cid/'calibration_inputs/parent_random10.csv.gz'
            if not parent.exists():
                if existing.exists() and existing.with_suffix('.json').exists():shutil.copy2(existing,parent);shutil.copy2(existing.with_suffix('.json'),parent_meta)
                else:
                    update(cid,'fetching_uncut_quality_parent');headroom(cid)
                    from adaptive_parent import fetch
                    fetch(parent,report['center_ra_deg'],report['center_dec_deg'],report['field_radius_deg']+.05)
            bright=inp/'bright_parent.csv.gz'
            if not bright.exists():
                update(cid,'fetching_bright_census')
                from selection_model import parent_query,download_parent
                headroom(cid);download_parent(parent_query(report['center_ra_deg'],report['center_dec_deg'],report['field_radius_deg']+.05,0.,12.,1),bright,100000)
            extension=root/'model_repairs'/cid;log=root/(cid+'_model.log')
            if not (extension/'bright_repair.json').exists():
                update(cid,'building_cluster_calibration');command(code/'bright_calibration.py',['--repair',base,'--cluster',cid,'--bright',bright,'--out',extension],log,cid)
            if not (extension/'preflight.json').exists():
                update(cid,'testing_model_inputs');command(code/'bright_calibration.py',['--finalize','--out',extension],log,cid)
            package=root/'packages'/(cid+'_portable_v1')
            if not package.with_suffix('.zip').exists():
                update(cid,'exporting_model_package')
                from portable_cluster import export
                headroom(cid);export(extension,cid,package,analysis_kind='model_comparison')
            validation=root/'package_validation'/cid
            if not (validation/'validation.json').exists():
                update(cid,'testing_relocated_package')
                command(HERE/'validate_transfer_package.py',['--bundle',package,'--simulation',extension/'preflight/simulation','--out',validation],log,cid)
            update(cid,'conditional_model_package_ready',archive=str(package.with_suffix('.zip')),mass_loss_ready=False,
                conditional_model_run_ready=True,physical_rate_measurement_validated=False,
                scope='Small real-input checks passed; full-grid and physical rate validation remain. R/T posterior is not computed.')
        except Exception as exc:update(cid,'model_preparation_failed',error=str(exc),mass_loss_ready=False)
    dump(root/'model_preparation_finished.json',dict(clusters=states,finished_utc=datetime.now(timezone.utc).isoformat()))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--retry-failed',action='store_true',help='Retry prior failures only after reviewing/correcting their cause')
    p.add_argument('--ready-index',type=Path,help='Explicitly reuse verified completed packages from this same transfer cohort');a=p.parse_args()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    run(a.root,a.retry_failed,a.ready_index)
