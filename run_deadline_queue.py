"""Continuous, resumable cluster queue: log errors and proceed to the next cluster."""
import argparse
import ctypes
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from empirical_population import measured_constraint
from run_cluster_survey import SAMPLE

HERE=Path(__file__).resolve().parent
REQUIRED=('coupled_dynamics.json','selection/selection.json','calibration_inputs/isochrone.dat',
          'science_calibration/bundle/calibration_manifest.json')


def memory_snapshot():
    """System RAM only; this is not the queue's private allocation."""
    if sys.platform=='win32':
        class Status(ctypes.Structure):
            _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(k,ctypes.c_ulonglong) for k in
                ('total','available','page_total','page_available','virtual_total','virtual_available','extended')]
        s=Status();s.length=ctypes.sizeof(s)
        if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(s)):
            return dict(total_gib=s.total/2**30,available_gib=s.available/2**30,used_gib=(s.total-s.available)/2**30)
    return {}


def read_json(path,default=None):
    return json.loads(path.read_text(encoding='utf-8-sig')) if path.exists() else default


def write_json(path,data):
    tmp=path.with_name(path.name+f'.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(data,indent=2,ensure_ascii=False),encoding='utf-8')
    for attempt in range(5):
        try:tmp.replace(path);return
        except PermissionError:
            if attempt==4:raise
            time.sleep(.1*(attempt+1))


def completed_reports(dest):
    for label in ('empirical_direct_v1','empirical_mf_v1'):
        reports=read_json(dest/label/'transfer_sensitivity.json',{}).get('reports',[])
        if len(reports)==3 and all(Path(p).is_file() for p in reports):return reports
    return []


def run(root,workers=12,*,retry_failed=False,memory_budget_gib=24.,clusters=None):
    if not 1<=workers<=12 or memory_budget_gib<=0:raise ValueError('Use 1–12 workers and a positive RAM budget')
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True)
    previous=read_json(root/'deadline_queue_status.json',{})
    old_entries={e['cluster']:e for e in previous.get('entries',[])}
    root=Path(root).resolve();rules=json.loads((HERE/'configs/survey_populations_v1.json').read_text())['clusters']
    for cid,extra in read_json(HERE/'configs/population_supplements_v1.json',{'clusters':{}})['clusters'].items():
        rules[cid]=dict(rules.get(cid,{}),**extra)
    html=HERE/'outputs/empirical_mf_source/parameter.html'
    ids=[cid for _,cid in SAMPLE]
    if clusters is not None:
        if set(clusters)-set(ids):raise ValueError('Unknown cluster in requested queue')
        ids=[c for c in ids if c in clusters]
    active={'checking_inputs','preparing_inputs','running_direct_empirical'}
    ids.sort(key=lambda c:(old_entries.get(c,{}).get('stage') not in active,not all((root/c/p).exists() for p in REQUIRED),
                           ['NGC5466','NGC5904','NGC362','NGC6144'].index(c) if c in ['NGC5466','NGC5904','NGC362','NGC6144'] else 99))
    status=dict(schema_version=2,started_utc=datetime.now(timezone.utc).isoformat(),stage='running',pid=os.getpid(),
        workers=workers,response_workers=workers,memory_budget_gib=memory_budget_gib,cluster_order=ids,entries=[],
        memory_policy='Soft scheduling budget; OS headroom check; no artificial RAM allocation')
    def event(kind,**data):
        with (root/'execution_events.jsonl').open('a',encoding='utf-8') as f:
            f.write(json.dumps(dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),event=kind,
                pid=os.getpid(),workers=workers,**data),ensure_ascii=False)+'\n')
    def save():
        status['updated_utc']=datetime.now(timezone.utc).isoformat()
        status['processed']=sum(e['stage'] not in active for e in status['entries']);status['total']=len(ids)
        write_json(root/'deadline_queue_status.json',status)
        lines=['# 連續批次執行紀錄','',f'更新 UTC：{status["updated_utc"]}',
            f'工作程序：{workers}；RAM 軟預算：{memory_budget_gib:g} GiB。',
            '完成或出錯後自動接下一個；失敗不視為科學結果，重啟接續檢查點。','',
            '| 星團 | 狀態 | 耗時（分鐘） | 原因或輸出 |','|---|---|---:|---|']
        for e in status['entries']:
            detail=str(e.get('reason',e.get('report',''))).replace('|','/').replace('\n',' ')
            lines.append(f"| {e['cluster']} | {e['stage']} | {e.get('elapsed_seconds',0)/60:.1f} | {detail[-1400:]} |")
        (root/'EXECUTION_LOG.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    def update(entry):
        entry['updated_utc']=datetime.now(timezone.utc).isoformat()
        status['updated_utc']=entry['updated_utc']
        save();event(entry['stage'],cluster=entry['cluster'],reason=entry.get('reason',''))
        print(entry['cluster'],entry['stage'],entry.get('reason',''),flush=True)
    def command(script,args,log,entry):
        snap=memory_snapshot()
        if snap and snap['available_gib']<max(1.,snap['total_gib']-memory_budget_gib):
            raise RuntimeError('Insufficient RAM headroom; defer this cluster')
        argv=[sys.executable,'-u',str(HERE/script),*map(str,args)]
        event('command_started',cluster=entry['cluster'],command=argv,log=str(log),system_ram=snap)
        with log.open('a',encoding='utf-8') as f:
            f.write('\n'+datetime.now(timezone.utc).isoformat()+' START '+script+'\n');f.flush()
            proc=subprocess.Popen(argv,cwd=HERE,
                stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW if sys.platform=='win32' else 0)
            entry['child_pid']=proc.pid;save()
            while True:
                try:code=proc.wait(timeout=30);break
                except subprocess.TimeoutExpired:
                    sample=memory_snapshot();entry['last_system_ram']=sample
                    entry['peak_system_used_gib']=max(entry.get('peak_system_used_gib',0),sample.get('used_gib',0))
                    save();event('resource_sample',cluster=entry['cluster'],child_pid=proc.pid,system_ram=sample)
        entry.pop('child_pid',None);event('command_finished',cluster=entry['cluster'],script=script,exit_code=code)
        if code:
            tail='\n'.join(log.read_text(encoding='utf-8',errors='replace').splitlines()[-8:])
            raise RuntimeError(f'{script} exited {code}: {tail}; log={log}')
    event('queue_started',memory_budget_gib=memory_budget_gib,retry_failed=retry_failed,cluster_order=ids,system_ram=memory_snapshot());save()
    for cid in ids:
        entry=dict(cluster=cid,stage='checking_inputs',physical_rate_measurement_validated=False)
        status['entries'].append(entry);dest=root/cid;started=time.monotonic()
        try:
            dest.mkdir(exist_ok=True)
            reports=completed_reports(dest)
            if reports:
                entry.update(stage='skipped_completed',report=reports[0]);update(entry);continue
            if cid=='Pal5':
                entry.update(stage='skipped_preserved',reason='Pal 5 rerun cancelled by user; preserve existing results');update(entry);continue
            old=old_entries.get(cid,{});old_job=read_json(dest/'empirical_job.json',{})
            if not retry_failed and (old.get('stage')=='blocked_or_failed' or old_job.get('stage')=='empirical_stage_failed'):
                entry.update(stage='blocked_or_failed',reason=old.get('reason') or old_job.get('error') or 'Previous failure; retry explicitly after correction',
                    skipped_previous_failure=True);update(entry);continue
            update(entry)
            rule=rules.get(cid,{})
            if not rule or rule.get('hold_reason'):raise ValueError(rule.get('hold_reason','Missing literature population age'))
            measured_constraint(html,cid)
            if not (dest/'selection/selection.json').exists():raise ValueError('Missing observed Gaia selection')
            if not all((dest/p).exists() for p in REQUIRED):
                entry['stage']='preparing_inputs';update(entry)
                command('run_remaining_survey.py',['--root',root,'--clusters',cid,'--prepare-only',
                    '--prepare-workers','1','--workers',workers,'--response-workers',workers],dest/'deadline_preparation.log',entry)
                state=json.loads((dest/'remaining_job.json').read_text())
                if state['stage']!='ready_conditional_grid':raise ValueError(state.get('error',state.get('reason',state['stage'])))
            entry['stage']='running_direct_empirical';update(entry)
            command('run_direct_empirical.py',['--root',root,'--cluster',cid,'--workers',workers],dest/'deadline_inference.log',entry)
            entry.update(stage='completed_conditional_posteriors_validation_pending',elapsed_seconds=time.monotonic()-started);update(entry)
        except Exception as exc:
            failure_stage=entry['stage']
            entry.update(stage='blocked_or_failed',failure_stage=failure_stage,reason=str(exc),elapsed_seconds=time.monotonic()-started);update(entry)
            try:
                job=read_json(dest/'empirical_job.json',{})
                if job:
                    job.update(stage='empirical_stage_failed',error=str(exc),updated_utc=entry['updated_utc'],workers=workers,
                        physical_rate_measurement_validated=False)
                    write_json(dest/'empirical_job.json',job)
            except (OSError,ValueError) as log_error:
                event('child_status_warning',cluster=cid,reason=str(log_error))
    status.update(stage='finished',finished_utc=datetime.now(timezone.utc).isoformat());save();event('queue_finished')
    return status


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=HERE/'outputs/cluster_survey_29_v1')
    p.add_argument('--workers',type=int,default=12)
    p.add_argument('--memory-budget-gib',type=float,default=24.)
    p.add_argument('--retry-failed',action='store_true');p.add_argument('--clusters',nargs='+');a=p.parse_args()
    if not 1<=a.workers<=12 or a.memory_budget_gib<=0:p.error('Use 1–12 workers and a positive RAM budget')
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
    run(a.root,a.workers,retry_failed=a.retry_failed,memory_budget_gib=a.memory_budget_gib,clusters=a.clusters)
