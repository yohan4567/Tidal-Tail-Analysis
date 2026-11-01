"""Isolated bright-support repair for NGC3201, with a serial optional local run."""
import argparse,json,os,shutil,subprocess,sys,time,zipfile
from pathlib import Path
from datetime import datetime,timezone
from provenance import sha256_file
from repair_four_populations import stage_code
from repair_parent_support import fetch_parent
from repair_seven_clusters import worker_lock
from transfer_status_io import write_json

HERE=Path(__file__).resolve().parent
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);write_json(p,d,status_file=p.name=='status.json')
def run(root,after=None):
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True);cid='NGC3201'
    def update(stage,**extra):
        d=dict(stage=stage,cluster=cid,updated_utc=datetime.now(timezone.utc).isoformat(),physical_rate_measurement_validated=False)
        d.update(extra);dump(root/'status.json',d);print(stage,extra.get('error',''),flush=True)
    def memory():
        from run_deadline_queue import memory_snapshot
        while (m:=memory_snapshot()) and m['available_gib']<8:update('waiting_for_RAM',available_gib=m['available_gib']);time.sleep(30)
    def command(script,args):
        memory()
        with (root/'prepare_model.log').open('a',encoding='utf-8') as log:
            r=subprocess.run([sys.executable,'-u',str(script),*map(str,args)],cwd=HERE,stdout=log,stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
        if r.returncode:raise RuntimeError('\n'.join((root/'prepare_model.log').read_text(encoding='utf-8',errors='replace').splitlines()[-10:]))
    try:
        from selection_bridge import DEFAULT_STUDY
        original=HERE/'outputs/cluster_survey_29_v1'/cid;study=root/'study';(study/'catalogs').mkdir(parents=True,exist_ok=True)
        for name in ('step3_tail_search.py','skymath.py','catalogs/master_clusters.csv'):
            if not (study/name).exists():shutil.copy2(DEFAULT_STUDY/name,study/name)
        source=root/'observed'/cid;source.mkdir(parents=True,exist_ok=True)
        if not (source/'selection').exists():shutil.copytree(original/'selection',source/'selection')
        for name in ('adopted_dynamics.json','coupled_dynamics.json'):
            if not (source/name).exists():shutil.copy2(original/name,source/name)
        inp=source/'calibration_inputs';inp.mkdir(exist_ok=True)
        for name in ('isochrone.dat','isochrone.json'):
            if not (inp/name).exists():shutil.copy2(original/'calibration_inputs'/name,inp/name)
        selection=read(source/'selection/selection.json')
        if selection['selection_file_sha256']!=sha256_file(study/'step3_tail_search.py'):raise ValueError('Selector differs from original')
        update('verifying_faint_parent')
        fetch_parent(inp/'parent_random10.csv.gz',selection,candidates=[original/'calibration_inputs/parent_random10.csv.gz'])
        update('fetching_bright_census');fetch_parent(inp/'bright_parent.csv.gz',selection,bright=True)
        os.environ['GC_TIDAL_STUDY']=str(study);base,code=stage_code(root,study)
        prepared=root/'prepared'
        if not (prepared/'bright_repair.json').exists():
            update('building_bright_calibration');command(code/'bright_calibration.py',
                ['--repair',base,'--cluster',cid,'--bright',inp/'bright_parent.csv.gz','--out',prepared])
        memory()
        from selection_model import QualityMap
        q=QualityMap.load(Path(read(prepared/'bright_repair.json')['reference']).parent/'quality')
        from recover_compact_repairs import stage
        checked=root/'checked';pooled=bool((q.arrays['n_available']<q.minimum).any())
        stage(prepared,checked,HERE/'outputs/rc2',pool=pooled)
        if not (checked/'repair_manifest.json').exists():
            update('testing_small_response',spatial_pooling=pooled);command(HERE/'recover_compact_repairs.py',['--finish',checked])
        manifest=read(checked/'repair_manifest.json');state=read(checked/'input_status.json')[cid]
        job=dict(cluster=cid,code=manifest['snapshot'],study=str(study),config=str(source/'coupled_dynamics.json'),
            selection=str(source/'selection'),isochrone=str(inp/'isochrone.dat'),reference=state['reference'],html=str(checked/'mf_source.html'),
            reuse_grid=str(original/'empirical_direct_v1/offset_0/grid'),cluster_age_myr=read(inp/'isochrone.json')['age_years']/1e6,
            frozen_sha256=manifest['frozen_sha256'],target_parameter_posterior_computed=False,physical_rate_measurement_validated=False)
        jobpath=root/'job.json'
        if jobpath.exists() and read(jobpath)!=job:raise ValueError('Job identity changed')
        dump(jobpath,job)
        update('testing_full_duration_support');command(HERE/'run_priority_model.py',
            ['--job',jobpath,'--out',root/'support','--workers',1,'--validate-only'])
        package=root/'pkgs/NGC3201_v2';validation=root/'v'
        if not package.with_suffix('.zip').exists():
            update('exporting_package');command(HERE/'portable_cluster.py',['export','--repair',checked,'--cluster',cid,
                '--out',package,'--analysis-kind','model_comparison'])
        if not (validation/'validation.json').exists():
            update('checking_relocated_response');command(HERE/'validate_transfer_package.py',
                ['--bundle',package,'--simulation',checked/'preflight/simulation','--out',validation])
        if read(validation/'validation.json')['stage']!='relocated_small_response_passed':raise ValueError('Relocation not validated')
        archive=package.with_suffix('.zip')
        with zipfile.ZipFile(archive) as z:
            if z.testzip() is not None:raise ValueError('ZIP integrity failed')
        entry=dict(archive=str(archive),bytes=archive.stat().st_size,archive_sha256=sha256_file(archive),zip_crc_verified=True,
            validation=str(validation/'validation.json'),analysis_kind='model_comparison',job=str(jobpath),physical_rate_measurement_validated=False)
        dump(root/'READY_PACKAGES.json',dict(ready_packages={cid:entry}));update('conditional_model_package_ready',**entry)
        if after is not None:
            terminal={'completed_model_comparisons_validation_pending','repaired_analysis_failed','inputs_not_ready'}
            while True:
                p=Path(after)/'local_dispatch_status.json';s=read(p).get('clusters',{}) if p.exists() else {}
                if all(s.get(c,{}).get('stage') in terminal for c in ('NGC5904','NGC6093')):break
                update('ready_waiting_for_M5_M80',**entry);time.sleep(30)
            update('running_repaired_model_comparisons',**entry)
            command(HERE/'run_priority_model.py',['--job',jobpath,'--out',root/'results','--workers',12])
            update('completed_model_comparisons_validation_pending',**entry)
    except Exception as exc:update('repair_failed',error=str(exc))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--run-after',type=Path);a=p.parse_args()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    with worker_lock(a.root):run(a.root,a.run_after)
