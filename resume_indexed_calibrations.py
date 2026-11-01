"""Resume failed dense-field parents with smaller indexed patches and isolated outputs."""
import argparse,json,os,shutil,subprocess,sys,time,zipfile
from pathlib import Path
from datetime import datetime,timezone
from provenance import sha256_file
from repair_parent_support import fetch_parent
from repair_seven_clusters import worker_lock
from transfer_status_io import write_json

HERE=Path(__file__).resolve().parent
SOURCES={'NGC6266':'seven_cluster_repairs_20261101_v1/calibration',
         'NGC6366':'other_pc_ten_20261101_v1','NGC6304':'other_pc_ten_20261101_v1'}
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);write_json(p,d,status_file=p.name=='status.json')
def run(root,clusters):
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True);states={}
    def update(cid,stage,**extra):
        values=dict(stage=stage,updated_utc=datetime.now(timezone.utc).isoformat(),physical_rate_measurement_validated=False)
        values.update(extra);states.setdefault(cid,{}).update(values)
        dump(root/'status.json',dict(clusters=states));print(cid,stage,extra.get('error',''),flush=True)
    def memory(cid):
        from run_deadline_queue import memory_snapshot
        while (m:=memory_snapshot()) and m['available_gib']<8:
            update(cid,'waiting_for_RAM',available_gib=m['available_gib']);time.sleep(30)
    def command(script,args,cid):
        memory(cid)
        with (root/(cid+'.log')).open('a',encoding='utf-8') as log:
            r=subprocess.run([sys.executable,'-u',str(script),*map(str,args)],cwd=HERE,stdout=log,stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
        if r.returncode:raise RuntimeError('\n'.join((root/(cid+'.log')).read_text(encoding='utf-8',errors='replace').splitlines()[-10:]))
    for cid in clusters:
        try:
            source_root=HERE/'outputs'/SOURCES[cid];base=source_root/'model_base';settings=read(base/'repair_manifest.json')
            source=Path(settings['source'])/cid;inp=source/'calibration_inputs';s=read(source/'selection/selection.json')
            os.environ['GC_TIDAL_STUDY']=settings['study']
            update(cid,'downloading_smaller_indexed_patches')
            faint=fetch_parent(inp/'parent_random10.csv.gz',s,prefer_tiles=True,tile_level=7,tiles_per_request=1,
                max_requests=1200,before_merge=lambda:memory(cid),reuse_dirs=[inp,inp/'parent_random10.csv.gz.parts'])
            update(cid,'downloading_bright_census',faint_parent_stars=faint['n_sources'])
            bright=inp/'bright_parent.csv.gz';fetch_parent(bright,s,bright=True)
            original=root/cid/'prepared'
            if not (original/'bright_repair.json').exists():
                update(cid,'building_calibration');command(Path(settings['snapshot'])/'bright_calibration.py',
                    ['--repair',base,'--cluster',cid,'--bright',bright,'--out',original],cid)
            memory(cid)
            from selection_model import QualityMap
            quality=QualityMap.load(Path(read(original/'bright_repair.json')['reference']).parent/'quality')
            pooled=bool((quality.arrays['n_available']<quality.minimum).any())
            from recover_compact_repairs import stage
            checked=root/cid/'checked'
            stage(original,checked,HERE/'outputs/rc2',pool=pooled)
            if not (checked/'repair_manifest.json').exists():
                update(cid,'testing_rebuilt_response',spatial_pooling=pooled)
                command(HERE/'recover_compact_repairs.py',['--finish',checked],cid)
            package=root/'pkgs'/(cid+'_v2');validation=root/'v'/cid
            if not package.with_suffix('.zip').exists():
                update(cid,'exporting_package');command(HERE/'portable_cluster.py',['export','--repair',checked,
                    '--cluster',cid,'--out',package,'--analysis-kind','model_comparison'],cid)
            if not (validation/'validation.json').exists():
                update(cid,'checking_relocated_response');command(HERE/'validate_transfer_package.py',
                    ['--bundle',package,'--simulation',checked/'preflight/simulation','--out',validation],cid)
            if read(validation/'validation.json')['stage']!='relocated_small_response_passed':raise ValueError('Relocation not validated')
            archive=package.with_suffix('.zip')
            with zipfile.ZipFile(archive) as z:
                if z.testzip() is not None:raise ValueError('ZIP integrity failed')
            entry=dict(archive=str(archive),bytes=archive.stat().st_size,archive_sha256=sha256_file(archive),zip_crc_verified=True,
                validation=str(validation/'validation.json'),analysis_kind='model_comparison',physical_rate_measurement_validated=False)
            globalpath=HERE/'outputs/other_pc_all_remaining_20261101_v1/READY_PACKAGES.json';index=read(globalpath)
            index['ready_packages'][cid]=entry;index['updated_utc']=datetime.now(timezone.utc).isoformat();dump(globalpath,index)
            local=read(root/'READY_PACKAGES.json') if (root/'READY_PACKAGES.json').exists() else dict(ready_packages={})
            local['ready_packages'][cid]=entry;dump(root/'READY_PACKAGES.json',local)
            update(cid,'conditional_model_package_ready',**entry)
        except Exception as exc:update(cid,'repair_failed',error=str(exc))
    dump(root/'finished.json',dict(clusters=states,finished_utc=datetime.now(timezone.utc).isoformat()))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--clusters',nargs='+',choices=list(SOURCES),default=list(SOURCES));a=p.parse_args()
    for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
    with worker_lock(a.root):run(a.root,a.clusters)
