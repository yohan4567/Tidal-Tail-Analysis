"""Run the six user-prioritized repairs serially, preserving original results."""
import argparse,json,os,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
from provenance import sha256_file

HERE=Path(__file__).resolve().parent
ORDER=['NGC6144','NGC5466','NGC5904','NGC6093','NGC2298','NGC3201']


def run(root,workers):
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True)
    names=('priority_queue.py','prepare_priority_job.py','run_priority_model.py','likelihood_comparison.py','bright_calibration.py')
    plan=dict(order=ORDER,workers=workers,source_sha256={n:sha256_file(HERE/n) for n in names},target_parameter_posterior_computed=False)
    p=root/'plan.json'
    if p.exists() and json.loads(p.read_text())!=plan:raise ValueError('Priority controller changed; choose a new version')
    p.write_text(json.dumps(plan,indent=2));states={}
    def update(cid,stage,**extra):
        states.setdefault(cid,{}).update(stage=stage,updated_utc=datetime.now(timezone.utc).isoformat(),**extra)
        (root/'status.json').write_text(json.dumps(dict(order=ORDER,workers=workers,entries=states,target_parameter_posterior_computed=False),indent=2))
        print(cid,stage,extra.get('error',''),flush=True)
    for cid in ORDER:update(cid,'queued')
    def command(script,args,log,cid):
        for name,sha in plan['source_sha256'].items():
            if sha256_file(HERE/name)!=sha:raise ValueError('Priority controller source changed: '+name)
        with log.open('a',encoding='utf-8') as f:
            f.write('\n'+datetime.now(timezone.utc).isoformat()+' '+script+'\n');f.flush()
            proc=subprocess.Popen([sys.executable,'-u',str(HERE/script),*map(str,args)],cwd=HERE,stdout=f,stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
            update(cid,states[cid]['stage'],child_pid=proc.pid)
            while True:
                try:code=proc.wait(timeout=30);break
                except subprocess.TimeoutExpired:
                    from run_deadline_queue import memory_snapshot
                    sample=memory_snapshot();update(cid,states[cid]['stage'],last_system_ram=sample)
                    with (root/'resource_samples.jsonl').open('a',encoding='utf-8') as r:r.write(json.dumps(dict(cluster=cid,updated_utc=datetime.now(timezone.utc).isoformat(),system_ram=sample))+'\n')
        if code:raise RuntimeError(script+': '+'\n'.join(log.read_text(encoding='utf-8',errors='replace').splitlines()[-8:]))
    for cid in ORDER:
        dest=root/cid;dest.mkdir(exist_ok=True)
        final=dest/'results/run_status.json'
        if final.exists() and json.loads(final.read_text()).get('stage')=='completed_model_comparisons_validation_pending':
            update(cid,'completed_model_comparisons_validation_pending');continue
        try:
            update(cid,'preparing_repair')
            command('prepare_priority_job.py',['--cluster',cid,'--out',dest/'prepared'],dest/'prepare.log',cid)
            update(cid,'validating_and_running_model_comparisons')
            command('run_priority_model.py',['--job',dest/'prepared/job.json','--out',dest/'results','--workers',workers],dest/'run.log',cid)
            update(cid,'completed_model_comparisons_validation_pending',physical_rate_measurement_validated=False)
        except Exception as exc:update(cid,'blocked_or_failed',error=str(exc),physical_rate_measurement_validated=False)
    (root/'finished.json').write_text(json.dumps(dict(finished_utc=datetime.now(timezone.utc).isoformat(),entries=states),indent=2))
    repair=HERE/'outputs/input_repairs_20261101_v1';override=repair/'dispatch_overrides.json'
    settings=json.loads(override.read_text()) if override.exists() else {'excluded_clusters':[]}
    excluded=set(settings.get('excluded_clusters',[]))
    if states.get('NGC3201',{}).get('stage')=='completed_model_comparisons_validation_pending':excluded.add('NGC3201')
    original=HERE/'outputs/cluster_survey_29_v1/NGC6171/empirical_direct_v1'
    if all((original/offset/'posterior/posterior_summary.json').exists() for offset in ('offset_0','offset_minus0p5','offset_minus1')):excluded.add('NGC6171')
    settings.update(excluded_clusters=sorted(excluded),updated_utc=datetime.now(timezone.utc).isoformat(),
        priority_completion_source=str(root/'finished.json'))
    override.write_text(json.dumps(settings,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--workers',type=int,default=12);a=p.parse_args()
    if not 1<=a.workers<=12:p.error('Use 1–12 workers locally')
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    run(a.root,a.workers)
