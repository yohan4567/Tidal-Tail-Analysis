"""Serial cluster queue with eight workers per cluster and explicit input blockers."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys

from empirical_population import measured_constraint
from run_cluster_survey import SAMPLE

HERE=Path(__file__).resolve().parent
REQUIRED=('coupled_dynamics.json','selection/selection.json','calibration_inputs/isochrone.dat',
          'science_calibration/bundle/calibration_manifest.json')


def run(root,workers=8):
    root=Path(root).resolve();rules=json.loads((HERE/'configs/survey_populations_v1.json').read_text())['clusters']
    html=HERE/'outputs/empirical_mf_source/parameter.html'
    # Current single-cluster launcher finishes NGC2298 before this queue starts.
    ids=[cid for _,cid in SAMPLE if cid not in ('Pal5','NGC288','NGC2298')]
    ids.sort(key=lambda c:(not all((root/c/p).exists() for p in REQUIRED),
                           ['NGC5466','NGC5904','NGC362','NGC6144'].index(c) if c in ['NGC5466','NGC5904','NGC362','NGC6144'] else 99))
    status=dict(started_utc=datetime.now(timezone.utc).isoformat(),workers=workers,entries=[])
    def update(entry):
        entry['updated_utc']=datetime.now(timezone.utc).isoformat()
        status['updated_utc']=entry['updated_utc']
        (root/'deadline_queue_status.json').write_text(json.dumps(status,indent=2))
        print(entry['cluster'],entry['stage'],entry.get('reason',''),flush=True)
    def command(script,args,log):
        with log.open('a',encoding='utf-8') as f:
            result=subprocess.run([sys.executable,'-u',str(HERE/script),*map(str,args)],cwd=HERE,
                stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW if sys.platform=='win32' else 0)
        if result.returncode:raise RuntimeError(f'{script} exited {result.returncode}; see {log}')
    for cid in ids:
        entry=dict(cluster=cid,stage='checking_inputs',physical_rate_measurement_validated=False)
        status['entries'].append(entry);update(entry);dest=root/cid;dest.mkdir(exist_ok=True)
        try:
            rule=rules.get(cid,{})
            if not rule or rule.get('hold_reason'):raise ValueError(rule.get('hold_reason','Missing literature population age'))
            measured_constraint(html,cid)
            if not (dest/'selection/selection.json').exists():raise ValueError('Missing observed Gaia selection')
            if not all((dest/p).exists() for p in REQUIRED):
                entry['stage']='preparing_inputs';update(entry)
                command('run_remaining_survey.py',['--root',root,'--clusters',cid,'--prepare-only',
                    '--prepare-workers','1','--workers',workers,'--response-workers',workers],dest/'deadline_preparation.log')
                state=json.loads((dest/'remaining_job.json').read_text())
                if state['stage']!='ready_conditional_grid':raise ValueError(state.get('error',state.get('reason',state['stage'])))
            entry['stage']='running_direct_empirical';update(entry)
            command('run_direct_empirical.py',['--root',root,'--cluster',cid,'--workers',workers],dest/'deadline_inference.log')
            entry['stage']='completed_conditional_posteriors_validation_pending';update(entry)
        except Exception as exc:
            entry.update(stage='blocked_or_failed',reason=str(exc));update(entry)
    status['finished_utc']=datetime.now(timezone.utc).isoformat()
    (root/'deadline_queue_status.json').write_text(json.dumps(status,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=HERE/'outputs/cluster_survey_29_v1')
    p.add_argument('--workers',type=int,default=8);a=p.parse_args()
    if not 1<=a.workers<=8:p.error('Use one to eight workers')
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
    run(a.root,a.workers)
