"""Adopt an exactly compatible existing dynamics configuration without loosening reuse."""
from pathlib import Path
import argparse,copy,json
from datetime import datetime,timezone
from provenance import sha256_file
from transfer_status_io import write_json

HERE=Path(__file__).resolve().parent
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def physical_config(config):
    value=copy.deepcopy(config)
    value.get('mass_evolution',{}).pop('source_description',None)
    return value
def prepare(root):
    root=Path(root).resolve();index=read(root/'READY_PACKAGES.json')
    for cid in ('NGC5904','NGC6093'):
        entry=index['ready_packages'][cid]
        if read(entry['validation'])['stage']!='relocated_small_response_passed' or not entry['zip_crc_verified']:
            raise ValueError('Validated repaired package required')
        oldjob=read(root/cid/'job.json');original=HERE/'outputs/cluster_survey_29_v1'/cid/'coupled_dynamics.json'
        requested=read(oldjob['config']);compatible=read(original)
        if physical_config(requested)!=physical_config(compatible):raise ValueError('Physical parameters differ; cannot reuse dynamics')
        job=copy.deepcopy(oldjob);job['config']=str(original)
        job['frozen_sha256'][str(original.resolve())]=sha256_file(original)
        job['exact_reuse_audit']=dict(original_configuration=str(original),original_sha256=sha256_file(original),
            repaired_configuration=oldjob['config'],repaired_sha256=sha256_file(oldjob['config']),
            only_difference='mass_evolution.source_description',all_physical_and_numerical_base_parameters_identical=True,
            policy='Retain the original descriptive annotation in a new local job; the strict exact-config reuse guard remains unchanged.')
        destination=root/cid/'local_v3';destination.mkdir(exist_ok=True)
        path=destination/'job.json'
        if path.exists() and read(path)!=job:raise ValueError('Existing retry job differs')
        write_json(path,job)
        entry['job']=str(path)
    write_json(root/'READY_PACKAGES.json',index)
    # The previous controller finished before this reconciliation. Only states
    # backed by real package validation are corrected; scientific artifacts stay fixed.
    if not (root/'finished.json').exists():raise ValueError('Preparation controller has not finished')
    status=read(root/'status.json')
    for cid,entry in index['ready_packages'].items():
        validation=read(entry['validation'])
        bundle=Path(entry['archive']).with_suffix('')
        if validation['stage']!='relocated_small_response_passed' or validation['bundle_manifest_sha256']!=sha256_file(bundle/'portable_manifest.json'):
            raise ValueError('Package validation identity differs')
        record=status['clusters'].setdefault(cid,{})
        record.pop('error',None);record.update(entry,stage='conditional_model_package_ready',
            updated_utc=datetime.now(timezone.utc).isoformat(),status_reconciled_after_controller_only_error=True)
    write_json(root/'status.json',status,status_file=True)
    write_json(root/'status_reconciliation.json',dict(reason='Duplicate status keyword after successful validation; ready registry, relocation identity and CRC flags checked.',
        clusters=list(index['ready_packages']),updated_utc=datetime.now(timezone.utc).isoformat()))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();prepare(a.root)
