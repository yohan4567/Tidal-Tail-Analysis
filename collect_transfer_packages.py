"""Collect only validated portable ZIPs as they become ready; preserve originals."""
import argparse,json,os,shutil,time,zipfile
from pathlib import Path
from datetime import datetime,timezone
from provenance import sha256_file
from repair_seven_clusters import worker_lock
from transfer_status_io import write_json


def read(path):
    try:return json.loads(Path(path).read_text(encoding='utf-8-sig'))
    except (FileNotFoundError,json.JSONDecodeError):return {}


def collect(entry,cid,out):
    original=Path(entry['archive']);bundle=original.with_suffix('')
    validation=read(entry['validation']);manifest=bundle/'portable_manifest.json'
    if (entry.get('analysis_kind')!='model_comparison' or not entry.get('zip_crc_verified')
            or validation.get('stage')!='relocated_small_response_passed'
            or validation.get('cluster')!=cid or validation.get('bundle_manifest_sha256')!=sha256_file(manifest)):
        raise ValueError('Missing or mismatched relocated model-comparison validation')
    digest=sha256_file(original)
    if digest!=entry['archive_sha256']:raise ValueError('Original archive changed')
    with zipfile.ZipFile(original) as z:
        if z.testzip() is not None:raise ValueError('Original ZIP CRC failed')
    destination=Path(out)/original.name
    if destination.exists():
        if sha256_file(destination)!=digest:raise ValueError('Existing transfer ZIP differs; choose a new transfer folder')
    else:
        temporary=destination.with_name(destination.name+f'.{os.getpid()}.partial')
        if temporary.exists():raise ValueError('An interrupted copy exists; inspect it before retry')
        shutil.copy2(original,temporary)
        if sha256_file(temporary)!=digest:raise ValueError('Transfer copy checksum failed')
        temporary.replace(destination)
    return dict(stage='validated_portable_package_copied',archive=str(destination.resolve()),bytes=destination.stat().st_size,
                sha256=digest,original_archive=str(original.resolve()),validation=str(Path(entry['validation']).resolve()),
                physical_rate_measurement_validated=False)


def run(owner,out,clusters):
    owner,out=Path(owner).resolve(),Path(out).resolve();states=read(out/'STATUS.json').get('clusters',{})
    for cid in clusters:states.setdefault(cid,dict(stage='waiting_for_validation'))
    while True:
        entries=read(owner/'READY_PACKAGES.json').get('ready_packages',{})
        source_states=read(owner/'status.json').get('clusters',{})
        finished=(owner/'finished.json').exists()
        for cid in clusters:
            if states[cid]['stage'] in ('validated_portable_package_copied','copy_failed','preparation_failed'):continue
            if cid in entries:
                try:states[cid]=collect(entries[cid],cid,out)
                except Exception as exc:states[cid]=dict(stage='copy_failed',error=str(exc))
            elif source_states.get(cid,{}).get('stage')=='repair_failed' or finished:
                states[cid]=dict(stage='preparation_failed',error=source_states.get(cid,{}).get('error','No validated package produced'))
        ready=sum(s['stage']=='validated_portable_package_copied' for s in states.values())
        terminal=all(s['stage'] in ('validated_portable_package_copied','copy_failed','preparation_failed') for s in states.values())
        stage='all_three_packages_ready' if ready==len(clusters) else ('finished_with_unavailable_packages' if terminal else 'collecting_validated_packages')
        write_json(out/'STATUS.json',dict(stage=stage,clusters=states,ready=ready,total=len(clusters),
            updated_utc=datetime.now(timezone.utc).isoformat(),physical_rate_measurement_validated=False),status_file=True)
        lines=['# 三顆星團移機包','',f'已通過檢查並複製：{ready}/{len(clusters)}。原始資料包保留原位。','',
            '| 星團 | 狀態 | ZIP |','|---|---|---|']
        for cid in clusters:
            s=states[cid];lines.append(f"| {cid} | {s['stage']} | {Path(s['archive']).name if 'archive' in s else '尚未可用'} |")
        lines += ['','## 在另一台 Windows 電腦執行','',
            '僅使用 STATUS.json 中標為 validated_portable_package_copied 的 ZIP。解壓縮後，在包含 portable_cluster.py 的資料夾開啟 PowerShell。',
            '需安裝 64 位元 Python 3.11。每顆資料包可獨立執行；建議逐顆運行，不要同時啟動三個完整分析。','',
            '```powershell','py -3.11 -m venv .venv','.\\.venv\\Scripts\\python.exe -m pip install -r requirements.txt',
            '.\\.venv\\Scripts\\python.exe portable_cluster.py verify',
            '.\\.venv\\Scripts\\python.exe -u portable_cluster.py run --workers 24 --response-workers 8',
            '.\\.venv\\Scripts\\python.exe portable_cluster.py status','```','',
            '24 workers 適用於先前提供的 24 核、64 GB RAM 電腦；觀測工作限制為 8。執行後保持資料夾位置不變。',
            '結果位於 results/offset_*/comparison/，包括 likelihood_grid.png 與 comparison.json。',
            '這是可開始條件式模型分析的資料包，並非已驗證的絕對質量流失率結果。']
        (out/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
        if terminal:return
        time.sleep(30)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--owner',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--clusters',nargs='+',default=['NGC6304','NGC6355','NGC6541']);a=p.parse_args()
    with worker_lock(a.out):run(a.owner,a.out,a.clusters)
