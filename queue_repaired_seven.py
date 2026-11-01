"""Append validated seven-cluster repairs after the existing local serial jobs."""
import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from provenance import sha256_file
from repair_seven_clusters import worker_lock
from transfer_status_io import write_json

HERE = Path(__file__).resolve().parent
TARGETS = ('NGC5904', 'NGC6093', 'NGC6341', 'NGC7099', 'NGC7078', 'NGC6139', 'NGC6266')
COMPLETE = 'completed_model_comparisons_validation_pending'
LOCAL_TERMINAL = {COMPLETE, 'repaired_analysis_failed', 'inputs_not_ready'}
TERMINAL = {COMPLETE, 'analysis_failed', 'inputs_failed', 'skipped_previous_failure'}


def read(path, default=None):
    try:
        return json.loads(Path(path).read_text(encoding='utf-8-sig'))
    except FileNotFoundError:
        return {} if default is None else default


def predecessors_finished(compact, bright):
    local = read(compact / 'local_dispatch_status.json').get('clusters', {})
    return (all(local.get(c, {}).get('stage') in LOCAL_TERMINAL for c in TARGETS[:2])
            and read(bright / 'status.json').get('stage') in {COMPLETE, 'repair_failed'})


def checked_entry(entry, cid):
    """Check package identity again before allowing an expensive local run."""
    archive = Path(entry['archive'])
    validation = read(entry['validation'])
    if (entry.get('analysis_kind') != 'model_comparison' or not entry.get('zip_crc_verified')
            or validation.get('stage') != 'relocated_small_response_passed'):
        raise ValueError('Package is not validated for model comparisons')
    if not archive.is_file() or sha256_file(archive) != entry['archive_sha256']:
        raise ValueError('Validated archive is missing or changed')
    bundle = archive.with_suffix('')
    if sha256_file(bundle / 'portable_manifest.json') != validation['bundle_manifest_sha256']:
        raise ValueError('Relocation validation does not match the package manifest')
    manifest = read(bundle / 'portable_manifest.json')
    if manifest.get('cluster') != cid or manifest.get('analysis_kind') != 'model_comparison':
        raise ValueError('Package cluster/mode differs')
    return bundle


def make_job(checked, cid, path):
    """Build an isolated local job from a newly validated network repair."""
    settings = read(checked / 'bright_repair.json')
    manifest = read(checked / 'repair_manifest.json')
    state = read(checked / 'input_status.json')[cid]
    source = Path(settings['source']) / cid
    inp = Path(settings['calibration_inputs_root'])
    job = dict(cluster=cid, code=settings['snapshot'], study=settings['study'],
               config=str(source / 'coupled_dynamics.json'), selection=str(source / 'selection'),
               isochrone=str(inp / 'isochrone.dat'), reference=state['reference'],
               html=str(checked / 'mf_source.html'), reuse_grid=None,
               cluster_age_myr=read(inp / 'isochrone.json')['age_years'] / 1e6,
               frozen_sha256=manifest['frozen_sha256'], target_parameter_posterior_computed=False,
               physical_rate_measurement_validated=False)
    if path.exists() and read(path) != job:
        raise ValueError('Existing queue job has changed')
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json(path, job)
    return path


def run(root, compact, network, bright, workers=12):
    root, compact, network, bright = (Path(p).resolve() for p in (root, compact, network, bright))
    plan = dict(clusters=list(TARGETS), workers=workers, compact_root=str(compact),
                network_root=str(network), wait_for_bright_retry=str(bright),
                analysis_kind='model_comparison', physical_rate_measurement_validated=False)
    if (root / 'plan.json').exists() and read(root / 'plan.json') != plan:
        raise ValueError('Queue plan differs; choose a new output root')
    write_json(root / 'plan.json', plan)
    states = read(root / 'status.json').get('clusters', {})
    for cid in TARGETS:
        states.setdefault(cid, dict(stage='queued_after_current_jobs'))

    def status(cid, stage, **extra):
        if cid is not None:
            states[cid].update(stage=stage, updated_utc=datetime.now(timezone.utc).isoformat(), **extra)
        write_json(root / 'status.json', dict(stage=stage if cid is None else 'processing_queue',
                   clusters=states, workers=workers, updated_utc=datetime.now(timezone.utc).isoformat(),
                   physical_rate_measurement_validated=False), status_file=True)
        if cid is not None:
            print(cid, stage, extra.get('error', ''), flush=True)

    while not predecessors_finished(compact, bright):
        status(None, 'waiting_for_existing_M5_M80_NGC3201')
        time.sleep(30)

    for cid in TARGETS:
        if states[cid]['stage'] in TERMINAL:
            continue
        # These two are owned by the already running dispatcher. Never duplicate it.
        if cid in TARGETS[:2]:
            previous = read(compact / 'local_dispatch_status.json')['clusters'][cid]
            if previous['stage'] == COMPLETE:
                final = Path(previous['result'])
                if read(final).get('stage') == COMPLETE:
                    status(cid, COMPLETE, result=str(final), reused_completed_run=True)
                    continue
            status(cid, 'skipped_previous_failure', error=previous.get('error', 'Previous run did not complete'))
            continue
        try:
            owner = network if cid == 'NGC6266' else compact
            while True:
                entry = read(owner / 'READY_PACKAGES.json').get('ready_packages', {}).get(cid)
                if entry is not None:
                    break
                preparation = read(owner / 'status.json').get('clusters', {}).get(cid, {})
                if preparation.get('stage') == 'repair_failed' or (owner / 'finished.json').exists():
                    status(cid, 'inputs_failed', error=preparation.get('error', 'No validated package produced'))
                    break
                status(cid, 'waiting_for_validated_inputs')
                time.sleep(30)
            if entry is None:
                continue
            checked_entry(entry, cid)
            job = (make_job(network / cid / 'checked', cid, root / cid / 'job.json')
                   if cid == 'NGC6266' else Path(entry['job']))
            out = root / cid / 'results'
            final = out / 'run_status.json'
            if read(final).get('stage') != COMPLETE:
                from run_deadline_queue import memory_snapshot
                while (m := memory_snapshot()) and m['available_gib'] < 8:
                    status(cid, 'waiting_for_RAM', available_gib=m['available_gib'])
                    time.sleep(30)
                status(cid, 'running_model_comparisons', job=str(job), result=str(final))
                log = root / (cid + '.log')
                with log.open('a', encoding='utf-8') as stream:
                    result = subprocess.run([sys.executable, '-u', str(HERE / 'run_priority_model.py'),
                        '--job', str(job), '--out', str(out), '--workers', str(workers)],
                        cwd=HERE, stdout=stream, stderr=subprocess.STDOUT,
                        creationflags=subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0)
                if result.returncode:
                    raise RuntimeError('Analysis failed; see ' + str(log))
            if read(final).get('stage') != COMPLETE:
                raise ValueError('Missing completed model-comparison marker')
            status(cid, COMPLETE, result=str(final))
        except Exception as exc:
            status(cid, 'analysis_failed', error=str(exc))
    status(None, 'queue_finished_check_individual_statuses')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--compact', type=Path, default=HERE / 'outputs/s7r2')
    parser.add_argument('--network', type=Path, default=HERE / 'outputs/netr2')
    parser.add_argument('--bright', type=Path, default=HERE / 'outputs/b3201')
    parser.add_argument('--workers', type=int, default=12)
    args = parser.parse_args()
    if not 1 <= args.workers <= 12:
        parser.error('Use 1–12 workers locally')
    for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[key] = '1'
    with worker_lock(args.root):
        run(args.root, args.compact, args.network, args.bright, args.workers)
