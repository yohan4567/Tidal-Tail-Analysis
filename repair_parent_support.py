"""Verified endpoint failover and reusable, uncut Gaia calibration parents."""
from pathlib import Path
from datetime import datetime,timezone
from urllib.error import HTTPError
import json,shutil,time
import numpy as np
from provenance import sha256_file
from selection_model import parent_query,download_parent

ENDPOINTS=('https://gaia.aip.de/tap/sync','https://gea.esac.esa.int/tap-server/tap/sync')

def matches(path,ra,dec,radius,gmin,gmax,modulus):
    path=Path(path);meta=path.with_suffix('.json')
    if not path.is_file() or not meta.is_file():return False
    m=json.loads(meta.read_text(encoding='utf-8-sig'))
    return (m.get('query_status')=='OK' and m.get('parent_is_unfiltered') is True
        and not m.get('synthetic_calibration',False) and m.get('g_min')==gmin and m.get('g_max')==gmax
        and m.get('random_sampling_modulus')==modulus and m.get('radius_deg')==radius
        and np.allclose(m.get('center',[np.inf,np.inf]),[ra,dec],rtol=0,atol=1e-8)
        and m.get('sha256')==sha256_file(path))

def fetch_parent(path,selection,*,bright=False,candidates=(),reuse_dirs=(),before_merge=None):
    from adaptive_parent import fetch,QueryTimeout
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    ra,dec=selection['center_ra_deg'],selection['center_dec_deg'];radius=selection['field_radius_deg']+.05
    gmin,gmax,modulus=(0.,12.,1) if bright else (12.,21.,10)
    query=parent_query(ra,dec,radius,gmin,gmax,modulus)
    if path.exists():
        if not matches(path,ra,dec,radius,gmin,gmax,modulus):raise ValueError('Invalid or incompatible parent checkpoint: '+str(path))
        return json.loads(path.with_suffix('.json').read_text())
    # The owning worker holds a run lock. Preserve an interrupted network body
    # as evidence rather than accepting or overwriting it as a complete parent.
    temporary=path.with_suffix(path.suffix+'.download.vot')
    if temporary.exists():temporary.rename(temporary.with_name(temporary.name+'.interrupted_'+str(time.time_ns())))
    for old in map(Path,candidates):
        if matches(old,ra,dec,radius,gmin,gmax,modulus):
            if bright and json.loads(old.with_suffix('.json').read_text())['query_adql']!=query:continue
            shutil.copy2(old,path);shutil.copy2(old.with_suffix('.json'),path.with_suffix('.json'))
            return json.loads(path.with_suffix('.json').read_text())
    attempts=path.parent/(path.name+'.attempts.jsonl')
    def download(q,out,limit,endpoint=None):
        order=ENDPOINTS if endpoint is None else (endpoint,*[ep for ep in ENDPOINTS if ep!=endpoint])
        failures=[];last=None
        for ep in order:
            for number in range(2):
                try:return download_parent(q,out,limit,ep)
                except (OSError,EOFError) as exc:
                    last=exc
                    entry=dict(utc=datetime.now(timezone.utc).isoformat(),endpoint=ep,attempt=number+1,
                               error=repr(exc),query_sha256=__import__('hashlib').sha256(q.encode()).hexdigest())
                    if isinstance(exc,HTTPError):
                        entry['response_excerpt']=exc.read(2400).decode('utf-8',errors='replace')
                    with attempts.open('a',encoding='utf-8') as f:f.write(json.dumps(entry)+'\n')
                    failures.append(entry)
                    if Path(out).exists() or Path(out).with_suffix('.json').exists():raise
                    if 'statement timeout' in entry.get('response_excerpt','').lower():
                        raise QueryTimeout('TAP server explicitly reported statement timeout; use smaller indexed patches') from exc
                    if number==0:time.sleep(1)
        raise RuntimeError('Gaia parent endpoints unavailable; see '+str(attempts)+'; '+str(failures[-1]['error'])) from last
    try:
        meta=download(query,path,100000)
        meta.update(center=[ra,dec],radius_deg=radius,g_min=gmin,g_max=gmax,random_sampling_modulus=modulus)
        path.with_suffix('.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    except RuntimeError as exc:
        if bright or ('OVERFLOW' not in str(exc) and not isinstance(exc,QueryTimeout)):raise
        meta=fetch(path,ra,dec,radius,g_min=gmin,g_max=gmax,modulus=modulus,
            endpoint=ENDPOINTS[0],reuse_dirs=reuse_dirs,downloader=download,before_merge=before_merge)
    if not matches(path,ra,dec,radius,gmin,gmax,modulus):raise ValueError('Final parent verification failed')
    return meta
