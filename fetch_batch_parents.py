"""Fetch cluster-specific unfiltered calibration parents; no population transfer."""
import json,subprocess,sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
root=Path(sys.argv[1]).resolve();here=Path(__file__).resolve().parent
entries=json.loads((root/'batch_audit.json').read_text())['entries']

def fetch(entry):
    cid=entry['cluster'];dest=root/cid;out=dest/'parent_random10.csv.gz';state=dest/'parent_download_status.json'
    if cid=='Pal5':return
    if out.exists() and out.with_suffix('.json').exists():return
    selection=json.loads((dest/'selection/selection.json').read_text())
    status=dict(cluster=cid,status='running',purpose='Quality/error calibration only; not member-selected input')
    state.write_text(json.dumps(status,indent=2))
    with (dest/'parent_download.log').open('w',encoding='utf-8') as log:
        for endpoint in ['esa','aip']:
            args=[sys.executable,'-u',str(here/'fetch_parent_tiles.py'),'--ra',str(selection['center_ra_deg']),
                '--dec',str(selection['center_dec_deg']),'--radius',str(selection['field_radius_deg']+.05),
                '--out',str(out),'--endpoint',endpoint]
            result=subprocess.run(args,cwd=here,stdout=log,stderr=subprocess.STDOUT)
            if result.returncode==0:break
    status.update(status='complete' if result.returncode==0 else 'failed',returncode=result.returncode)
    state.write_text(json.dumps(status,indent=2));print(cid,status['status'],flush=True)
with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(fetch,entries))
