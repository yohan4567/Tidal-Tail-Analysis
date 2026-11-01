"""Prioritize cached targets without waiting for unrelated slow field downloads."""
import argparse,json,os,shutil,time
from pathlib import Path
from datetime import datetime,timezone
from provenance import sha256_file

HERE=Path(__file__).resolve().parent

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,d):p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2),encoding='utf-8')

def run(root,source,clusters,base_snapshot=None,bright_parent=None):
    root,source=Path(root).resolve(),Path(source).resolve();root.mkdir(parents=True,exist_ok=True)
    info=read(source/'status.json')['clusters'];old=HERE/'outputs/input_repairs_20261101_v1'
    # Only admit targets with a complete observed selection and verified cached
    # faint parents. No extra sky survey is placed in front of these packages.
    for cid in clusters:
        parent=old/cid/'calibration_inputs/parent_random10.csv.gz'
        if not info.get(cid,{}).get('observation_ready'):raise ValueError('Observed selection not ready: '+cid)
        meta=read(parent.with_suffix('.json'))
        if meta.get('query_status')!='OK' or meta['sha256']!=sha256_file(parent):raise ValueError('Unverified cached parent: '+cid)
    study=root/'study'
    for name in ('step3_tail_search.py','skymath.py','catalogs/master_clusters.csv'):
        dest=study/name;dest.parent.mkdir(parents=True,exist_ok=True)
        if not dest.exists():shutil.copy2(source/'study'/name,dest)
    fields=study/'gaia_fields';fields.mkdir(exist_ok=True)
    for cid in clusters:
        for suffix in ('_field.csv.gz','_field.json','_missing_tiles.json'):
            p=source/'study/gaia_fields'/(cid+suffix)
            if p.exists() and not (fields/p.name).exists():shutil.copy2(p,fields/p.name)
        dest=root/'observed'/cid
        if not dest.exists():
            dest.mkdir(parents=True)
            shutil.copytree(source/'observed'/cid/'selection',dest/'selection')
            shutil.copy2(source/'observed'/cid/'adopted_dynamics.json',dest/'adopted_dynamics.json')
    if bright_parent:
        if len(clusters)!=1:raise ValueError('A cached bright parent requires one target')
        bright_parent=Path(bright_parent).resolve();meta=read(bright_parent.with_suffix('.json'))
        if meta.get('query_status')!='OK' or meta['sha256']!=sha256_file(bright_parent):raise ValueError('Unverified bright parent')
        dest=root/'observed'/clusters[0]/'calibration_inputs';dest.mkdir(exist_ok=True)
        shutil.copy2(bright_parent,dest/'bright_parent.csv.gz');shutil.copy2(bright_parent.with_suffix('.json'),dest/'bright_parent.csv.json')
    if base_snapshot and not (root/'model_base/code').exists():
        original=read(Path(base_snapshot)/'repair_manifest.json')
        for name,sha in original['frozen_sha256'].items():
            if sha256_file(name)!=sha:raise ValueError('Source base changed: '+name)
        base=root/'model_base';code=base/'code'
        shutil.copytree(original['snapshot'],code,ignore=shutil.ignore_patterns('__pycache__','outputs'))
        # Only lossless cache compression changes. CSV values and all population,
        # selection, numerical-integration and calibration settings are retained.
        for name in ('bright_calibration.py','prepare_science_calibration.py'):
            p=code/name;text=p.read_text(encoding='utf-8')
            if "compression='gzip'" not in text:raise ValueError('Expected CSV compression call: '+name)
            p.write_text(text.replace("compression='gzip'","compression={'method':'gzip','compresslevel':1}"),encoding='utf-8')
        shutil.copy2(Path(base_snapshot)/'mf_source.html',base/'mf_source.html')
        files=[*code.glob('*.py'),*list((code/'configs').rglob('*.json')),base/'mf_source.html']
        dump(base/'repair_manifest.json',dict(source=str(root/'observed'),snapshot=str(code),study=str(study),
            source_manifest_sha256=sha256_file(Path(base_snapshot)/'repair_manifest.json'),
            cache_compression='lossless_gzip_level_1',preparation_source_sha256=sha256_file(Path(__file__)),
            frozen_sha256={str(p.resolve()):sha256_file(p) for p in files}))
    plan=dict(clusters=clusters,analysis_kind='model_comparison',priority='cached_inputs_first',
        source_observations=str(source),source_sha256=sha256_file(Path(__file__)),
        base_snapshot=str(Path(base_snapshot).resolve()) if base_snapshot else None,
        bright_parent_sha256=sha256_file(bright_parent) if bright_parent else None,
        bright_query_primary_endpoint='https://gaia.aip.de/tap/sync',physical_rate_measurement_validated=False)
    p=root/'transfer_plan.json'
    if p.exists() and read(p)!=plan:raise ValueError('Fast transfer plan changed; choose a new output version')
    dump(p,plan);dump(root/'status.json',dict(clusters={c:info[c] for c in clusters}))
    os.environ['GC_TIDAL_STUDY']=str(study)
    import selection_model
    original=selection_model.download_parent
    def download(query,out,max_records=250000,endpoint=None):
        # Same ADQL and completeness/overflow guards; endpoint is recorded by
        # the downloader. Failed attempts cannot leave an accepted truncated file.
        out=Path(out);last=None
        for ep in ('https://gaia.aip.de/tap/sync','https://gea.esac.esa.int/tap-server/tap/sync'):
            try:return original(query,out,max_records,ep)
            except (OSError,TimeoutError) as exc:
                last=exc
                if out.exists() or out.with_suffix('.json').exists():raise
                print('Bright query unavailable:',ep,str(exc),flush=True)
        raise RuntimeError('Both bright-query endpoints failed; no partial parent accepted') from last
    selection_model.download_parent=download
    from prepare_transfer_models import run as prepare
    prepare(root)
    result=read(root/'model_preparation_finished.json')
    result.update(finished_utc=datetime.now(timezone.utc).isoformat(),priority='cached_inputs_first')
    dump(root/'FAST_TRANSFER_RESULT.json',result)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--clusters',nargs='+',default=['NGC6397','NGC6121'])
    p.add_argument('--base-snapshot',type=Path,help='Verified transfer base; use faster lossless CSV cache compression in a new snapshot')
    p.add_argument('--bright-parent',type=Path,help='Verified completed bright census for a single target');a=p.parse_args()
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
    run(a.root,a.source,a.clusters,a.base_snapshot,a.bright_parent)
