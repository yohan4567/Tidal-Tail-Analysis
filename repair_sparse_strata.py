"""Repair NGC6355 with an uncut G=12..19 census in an isolated version."""
import argparse
import json
import os
import subprocess
import sys
import time
import traceback
import zipfile
from datetime import datetime, timezone
from pathlib import Path

from provenance import sha256_file
from repair_seven_clusters import worker_lock
from transfer_status_io import write_json

HERE = Path(__file__).resolve().parent
CID = 'NGC6355'


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def run(root):
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    state = {}

    def status(stage, **extra):
        state.update(stage=stage, updated_utc=datetime.now(timezone.utc).isoformat(),
                     physical_rate_measurement_validated=False, **extra)
        if stage != 'repair_failed':
            state.pop('error', None)
        write_json(root/'status.json', dict(clusters={CID: state}), status_file=True)
        print(stage, extra, flush=True)

    def memory():
        from run_deadline_queue import memory_snapshot
        previous = state.get('stage', 'preparing')
        waited = False
        while (m := memory_snapshot()) and m['available_gib'] < 8:
            waited = True
            status('waiting_for_RAM', available_gib=m['available_gib'], required_gib=8)
            time.sleep(30)
        if waited:
            status(previous, available_gib=m['available_gib'] if m else None)

    def command(script, args):
        memory()
        with (root/'preparation.log').open('a', encoding='utf-8') as log:
            result = subprocess.run([sys.executable, '-u', str(script), *map(str,args)],
                cwd=HERE, stdout=log, stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
        if result.returncode:
            raise RuntimeError('\n'.join((root/'preparation.log').read_text(
                encoding='utf-8', errors='replace').splitlines()[-16:]))

    try:
        base = HERE/'outputs/other_pc_ten_20261101_v1/model_base'
        settings = read(base/'repair_manifest.json')
        source = Path(settings['source'])/CID
        selection = read(source/'selection/selection.json')
        os.environ['GC_TIDAL_STUDY'] = settings['study']
        supplement = root/'inputs/G12_19_census.csv.gz'
        status('downloading_uncut_G12_19_census')
        from repair_parent_support import fetch_parent
        fetch_parent(supplement, selection, magnitude_range=(12.,19.), sampling_modulus=1,
            prefer_tiles=True, tile_level=7, tiles_per_request=1, max_requests=3000,
            before_merge=memory)
        prepared = root/'prepared'
        if not (prepared/'bright_repair.json').exists():
            status('rebuilding_calibration_with_census')
            command(HERE/'bright_calibration.py', ['--repair',base,'--cluster',CID,
                '--bright',source/'calibration_inputs/bright_parent.csv.gz',
                '--supplement',supplement,'--out',prepared])
        memory()
        status('checking_training_support_and_heldout_quality')
        from recover_compact_repairs import stage
        checked = root/'checked'
        stage(prepared, checked, HERE/'outputs/rc4', pool=True,
              pool_screening='exact_beta_binomial_v2')
        if not (checked/'repair_manifest.json').exists():
            status('testing_real_input_response')
            command(HERE/'recover_compact_repairs.py', ['--finish',checked])
        package = root/'pkgs'/f'{CID}_v3'
        validation = root/'v'/CID
        if not package.with_suffix('.zip').exists():
            status('exporting_portable_package')
            command(HERE/'portable_cluster.py', ['export','--repair',checked,'--cluster',CID,
                '--out',package,'--analysis-kind','model_comparison'])
        if not (validation/'validation.json').exists():
            status('testing_relocated_response')
            command(HERE/'validate_transfer_package.py', ['--bundle',package,
                '--simulation',checked/'preflight/simulation','--out',validation])
        result = read(validation/'validation.json')
        if result.get('stage')!='relocated_small_response_passed' or result.get('cluster')!=CID:
            raise ValueError('Missing relocated NGC6355 response validation')
        archive = package.with_suffix('.zip')
        with zipfile.ZipFile(archive) as z:
            if z.testzip() is not None:
                raise ValueError('ZIP CRC verification failed')
        entry = dict(archive=str(archive), archive_sha256=sha256_file(archive),
            bytes=archive.stat().st_size, zip_crc_verified=True,
            validation=str(validation/'validation.json'), analysis_kind='model_comparison',
            physical_rate_measurement_validated=False)
        # Use a private transfer folder: the older collector still owns its index.
        from collect_transfer_packages import collect
        transfer = root/'transfer'
        transfer.mkdir(exist_ok=True)
        copied = collect(entry, CID, transfer)
        write_json(transfer/'STATUS.json', dict(clusters={CID:copied},ready=1,total=1))
        write_json(root/'READY_PACKAGES.json', dict(ready_packages={CID:entry}))
        indexpath = HERE/'outputs/other_pc_all_remaining_20261101_v1/READY_PACKAGES.json'
        index = read(indexpath)
        index['ready_packages'][CID] = entry
        index['updated_utc'] = datetime.now(timezone.utc).isoformat()
        write_json(indexpath,index)
        status('conditional_model_package_ready', **entry, transfer_archive=copied['archive'])
    except Exception as exc:
        (root/'failure_traceback.log').write_text(traceback.format_exc(),encoding='utf-8')
        status('repair_failed', error=str(exc))
        return 1
    return 0


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
        os.environ[key]='1'
    with worker_lock(args.root):
        sys.exit(run(args.root))
