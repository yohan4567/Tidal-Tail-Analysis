"""Complete Gaia parents with recursive overflow splitting and verified cache reuse."""
from http.client import IncompleteRead
from pathlib import Path
from datetime import datetime,timezone
from urllib.error import HTTPError,URLError
import hashlib,json,os,time
import numpy as np
import pandas as pd
from fetch_parent_tiles import tile_query
from selection_model import download_parent
from selection_bridge import DEFAULT_STUDY,load_study
from provenance import sha256_file

ENDPOINTS=('https://gea.esac.esa.int/tap-server/tap/sync','https://gaia.aip.de/tap/sync')


def split_patch(pixels,level):
    pixels=np.asarray(pixels,dtype=np.int64)
    if len(pixels)>1:
        mid=len(pixels)//2;return [(pixels[:mid],level),(pixels[mid:],level)]
    if level>=12:raise RuntimeError('Overflow persists at the finest supported sky tile')
    return [(np.array([4*int(pixels[0])+k],dtype=np.int64),level+1) for k in range(4)]


def fetch(out,ra,dec,radius,*,g_min=12.,g_max=21.,modulus=10,level=6,tiles_per_request=4,
          endpoint=ENDPOINTS[0],reuse_dirs=None,pixels=None,downloader=None,max_requests=3000):
    if not (0<=ra<360 and -90<=dec<=90 and 0<radius<=10 and g_min<g_max and 1<=level<=12 and tiles_per_request>=1 and modulus>=1):
        raise ValueError('Invalid parent geometry, sampling or tiling')
    out=Path(out);metadata=out.with_suffix('.json')
    if out.exists() or metadata.exists():raise ValueError('Final parent output already exists')
    parts=out.parent/(out.name+'.parts');parts.mkdir(parents=True,exist_ok=True)
    downloader=downloader or download_parent
    if pixels is None:
        T=load_study(DEFAULT_STUDY,True);pixels=np.sort(T.skymath.cone_tiles(level,ra,dec,radius))
    pixels=np.asarray(pixels,dtype=np.int64)
    if not len(pixels) or len(set(pixels.tolist()))!=len(pixels):raise ValueError('Nonempty disjoint sky tiles required')
    reuse_dirs=reuse_dirs or json.loads(os.environ.get('GC_PARENT_REUSE_DIRS','[]'))
    cache={};leaves=[];requests=0;progress=out.parent/'adaptive_parent_status.json'
    for directory in [parts,*map(Path,reuse_dirs)]:
        if not directory.exists():continue
        for meta in directory.rglob('*.csv.json'):
            data=json.loads(meta.read_text(encoding='utf-8'));file=meta.with_suffix('.gz')
            if file.exists() and data.get('query_status')=='OK' and data.get('query_adql'):
                cache[data['query_adql']]=(file,meta,data)
    def save(stage,**extra):
        progress.write_text(json.dumps(dict(stage=stage,updated_utc=datetime.now(timezone.utc).isoformat(),
            requests=requests,completed_leaf_patches=len(leaves),**extra),indent=2),encoding='utf-8')
    def leaf(group,order):
        nonlocal requests
        query=tile_query(group,order,g_min,g_max,modulus)
        if query in cache:
            file,meta,data=cache[query]
            if sha256_file(file)!=data['sha256']:raise ValueError('Cached parent checksum mismatch: '+str(file))
            leaves.append((file,meta,group,order,query,True));save('downloading',cache_hit=str(file));return
        stem=hashlib.sha256(query.encode()).hexdigest()[:24];file=parts/(stem+'.csv.gz');meta=file.with_suffix('.json')
        partial=file.with_suffix(file.suffix+'.download.vot')
        if partial.exists():partial.rename(partial.with_name(partial.name+'.interrupted_'+str(time.time_ns())))
        selected=[endpoint,*[ep for ep in ENDPOINTS if ep!=endpoint]];last=None
        for ep in selected:
            for attempt in range(2):
                if requests>=max_requests:raise RuntimeError('Parent request budget exceeded; resume after investigating the log')
                requests+=1;save('downloading',level=order,pixels=group.tolist(),endpoint=ep)
                try:
                    downloader(query,file,100000,ep)
                    data=json.loads(meta.read_text(encoding='utf-8'))
                    if data.get('query_status')!='OK' or sha256_file(file)!=data['sha256']:raise ValueError('Unverified parent patch')
                    leaves.append((file,meta,group,order,query,False));save('downloading',last_patch=str(file));return
                except RuntimeError as exc:
                    if 'OVERFLOW' not in str(exc):raise
                    save('splitting_overflow',level=order,pixels=group.tolist())
                    for child,child_order in split_patch(group,order):leaf(child,child_order)
                    return
                except (TimeoutError,OSError,EOFError,IncompleteRead,URLError) as exc:
                    if isinstance(exc,HTTPError) and exc.code not in (408,429,500,502,503,504):raise
                    last=exc
                    if attempt==0:time.sleep(2)
        raise RuntimeError('Both Gaia endpoints failed for a patch; no truncated parent accepted') from last
    try:
        for start in range(0,len(pixels),tiles_per_request):leaf(pixels[start:start+tiles_per_request],level)
        save('merging');frames=[];patches=[]
        for file,meta,group,order,query,reused in leaves:
            data=json.loads(meta.read_text(encoding='utf-8'));frame=pd.read_csv(file,low_memory=False)
            if 'source_id' not in frame and 'datalinkid' in frame:frame=frame.rename(columns={'datalinkid':'source_id'})
            fac=2**35*4**(12-order)
            if 'source_id' not in frame or not pd.api.types.is_integer_dtype(frame.source_id) or not (frame.source_id//fac).isin(group).all():
                raise ValueError('Inexact or out-of-tile source IDs: '+str(file))
            frames.append(frame);patches.append(dict(path=str(file.resolve()),metadata=str(meta.resolve()),sha256=data['sha256'],
                level=order,pixels=group.tolist(),query_adql=query,endpoint=data.get('endpoint'),reused=reused,n_sources=len(frame)))
        merged=pd.concat(frames,ignore_index=True)
        if merged.source_id.duplicated().any():raise ValueError('Duplicate sources across disjoint leaf patches')
        a,d,a0,d0=map(np.radians,(merged.ra,merged.dec,ra,dec));sep=np.degrees(np.arccos(np.clip(np.sin(d0)*np.sin(d)+np.cos(d0)*np.cos(d)*np.cos(a-a0),-1,1)))
        merged=merged.loc[sep<=radius].copy();merged.to_csv(out,index=False,compression='infer')
        data=dict(parent_is_unfiltered=True,synthetic_calibration=False,query_status='OK',
            source_description='Gaia DR3 complete disjoint HEALPix leaf queries, recursively split only on overflow; sky/G and independent random_index sampling; exact cone retained locally; no quality/member cuts',
            query_adql=tile_query(pixels,level,g_min,g_max,modulus),patch_queries=[p['query_adql'] for p in patches],
            endpoint=endpoint,endpoints_used=sorted(set(p['endpoint'] for p in patches if p['endpoint'])),
            center=[ra,dec],radius_deg=radius,g_min=g_min,g_max=g_max,random_sampling_modulus=modulus,
            random_sampling_fraction=1/modulus,n_sources=len(merged),sha256=sha256_file(out),leaf_patches=patches,
            downloader_sha256=sha256_file(Path(__file__)))
        metadata.write_text(json.dumps(data,indent=2),encoding='utf-8');save('complete',n_sources=len(merged),output=str(out.resolve()))
        return data
    except Exception as exc:
        save('failed',error=str(exc));raise


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('ra','dec','radius'):p.add_argument('--'+key,type=float,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--sampling-modulus',type=int,default=10)
    p.add_argument('--tiles-per-request',type=int,default=4)
    p.add_argument('--endpoint',choices=('esa','aip'),default='esa')
    a=p.parse_args()
    result=fetch(a.out,a.ra,a.dec,a.radius,modulus=a.sampling_modulus,
        tiles_per_request=a.tiles_per_request,endpoint=ENDPOINTS[a.endpoint=='aip'])
    print(json.dumps(dict(output=str(a.out),n_sources=result['n_sources']),indent=2))
