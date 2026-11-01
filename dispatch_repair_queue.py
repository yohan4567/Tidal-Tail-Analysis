"""Apply scheduling exclusions independently of rewritten preparation status."""
import argparse,importlib.util,json,os,sys,subprocess
from pathlib import Path
from datetime import datetime,timezone
from provenance import sha256_file


def selected_clusters(available,overrides):
    excluded=overrides.get('excluded_clusters',[])
    if not isinstance(excluded,list) or any(not isinstance(c,str) for c in excluded) or set(excluded)-set(available):
        raise ValueError('Exclusions must name clusters in the frozen repair set')
    return tuple(c for c in available if c not in excluded)


def run(repair,workers):
    repair=Path(repair).resolve();manifest=repair/'repair_manifest.json'
    info=json.loads(manifest.read_text(encoding='utf-8'))
    for file,digest in info['frozen_sha256'].items():
        if sha256_file(file)!=digest:raise ValueError('Frozen repair input changed: '+file)
    override=repair/'dispatch_overrides.json'
    settings=json.loads(override.read_text(encoding='utf-8')) if override.exists() else {}
    if settings.get('ready_local_repairs_root'):
        # The existing scheduler calls this dispatcher after its priority pass.
        # An explicit registry adds validated retries without changing a live grid.
        subprocess.run([sys.executable,'-u',str(Path(__file__).resolve().parent/'dispatch_ready_local_repairs.py'),
            '--root',settings['ready_local_repairs_root'],'--workers',str(workers)],check=True)
    code=Path(info['snapshot']);sys.path.insert(0,str(code))
    spec=importlib.util.spec_from_file_location('_frozen_repair_dispatch',code/'repair_survey_inputs.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    selected=selected_clusters(module.CLUSTERS,settings)
    # Only the dispatch list changes. Frozen numerical/model files and their
    # identities remain unchanged; the scientific runner verifies them again.
    module.CLUSTERS=selected
    record=dict(selected_clusters=list(selected),excluded_clusters=settings.get('excluded_clusters',[]),
        source_snapshot=str(code),manifest_sha256=sha256_file(manifest),dispatch_sha256=sha256_file(Path(__file__)),
        overrides_sha256=sha256_file(override) if override.exists() else None,
        started_utc=datetime.now(timezone.utc).isoformat(),workers=workers)
    (repair/'dispatch_audit.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    module.run_repaired(repair,workers)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repair',type=Path,required=True);p.add_argument('--workers',type=int,default=12);a=p.parse_args()
    if not 1<=a.workers<=12:p.error('Use 1–12 workers for this local repair queue')
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    run(a.repair,a.workers)
