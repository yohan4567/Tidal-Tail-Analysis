"""Serially retry M5/M80 only after their repaired inputs pass validation."""
import argparse,json,os,subprocess,sys,time
from pathlib import Path
from datetime import datetime,timezone
from repair_seven_clusters import worker_lock
from transfer_status_io import write_json

HERE=Path(__file__).resolve().parent
TARGETS=('NGC5904','NGC6093')

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,v):write_json(Path(p),v,status_file=Path(p).name=='local_dispatch_status.json')
def preparing(root):
    try:
        with worker_lock(root):return False
    except (PermissionError,OSError):return True

def run(root,workers):
    root=Path(root).resolve();states={}
    def status(cid,stage,**extra):
        states.setdefault(cid,{}).update(stage=stage,updated_utc=datetime.now(timezone.utc).isoformat(),**extra)
        write(root/'local_dispatch_status.json',dict(clusters=states,workers=workers,
            physical_rate_measurement_validated=False));print(cid,stage,flush=True)
    for cid in TARGETS:
        while True:
            ready=read(root/'READY_PACKAGES.json').get('ready_packages',{}) if (root/'READY_PACKAGES.json').exists() else {}
            entry=ready.get(cid)
            if entry is not None:break
            failed=read(root/'status.json').get('clusters',{}).get(cid,{}) if (root/'status.json').exists() else {}
            if failed.get('stage')=='repair_failed' or not preparing(root):
                status(cid,'inputs_not_ready',error=failed.get('error','Preparation has not produced a validated ready package'))
                break
            status(cid,'waiting_for_validated_repair');time.sleep(30)
        if entry is None:continue
        try:
            if not entry.get('zip_crc_verified') or read(entry['validation'])['stage']!='relocated_small_response_passed':
                raise ValueError('Repaired package validation is incomplete')
            job=Path(entry['job']);out=job.parent/'results'
            final=out/'run_status.json'
            if not final.exists() or read(final).get('stage')!='completed_model_comparisons_validation_pending':
                status(cid,'running_repaired_model_comparisons',job=str(job))
                with (job.parent/'production.log').open('a',encoding='utf-8') as log:
                    result=subprocess.run([sys.executable,'-u',str(HERE/'run_priority_model.py'),
                        '--job',str(job),'--out',str(out),'--workers',str(workers)],cwd=HERE,stdout=log,stderr=subprocess.STDOUT,
                        creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
                if result.returncode:raise RuntimeError('Repaired analysis failed; see '+str(job.parent/'production.log'))
            if read(final)['stage']!='completed_model_comparisons_validation_pending':raise ValueError('Missing completed model-comparison marker')
            status(cid,'completed_model_comparisons_validation_pending',result=str(final))
            # The later general queue must not start the original broken version.
            override=HERE/'outputs/cluster_survey_29_v1/dispatch_overrides.json';settings=read(override)
            settings['excluded_clusters']=sorted(set(settings.get('excluded_clusters',[]))|{cid})
            settings.setdefault('completed_repaired_local_runs',{})[cid]=str(final)
            settings['updated_utc']=datetime.now(timezone.utc).isoformat();write(override,settings)
        except Exception as exc:status(cid,'repaired_analysis_failed',error=str(exc))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--workers',type=int,default=12);a=p.parse_args()
    if not 1<=a.workers<=12:p.error('Use 1–12 workers locally')
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    run(a.root,a.workers)
