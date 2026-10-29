"""Resumable unfiltered Gaia parent: indexed HEALPix source-id tiles, small requests."""
import argparse
from http.client import IncompleteRead
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from provenance import sha256_file
from selection_bridge import DEFAULT_STUDY,load_study
from selection_model import PARENT_COLUMNS,download_parent


def tile_query(pixels,level,g_min,g_max,modulus):
    fac=2**35*4**(12-level)
    ranges=' OR '.join(f'(source_id BETWEEN {int(p)*fac} AND {(int(p)+1)*fac-1})' for p in pixels)
    q=f'SELECT {PARENT_COLUMNS} FROM gaiadr3.gaia_source WHERE ({ranges}) AND phot_g_mean_mag >= {g_min} AND phot_g_mean_mag < {g_max}'
    if modulus>1: q+=f' AND MOD(random_index,{modulus}) = 0'
    return q


def fetch(out,ra,dec,radius,*,g_min=12.,g_max=21.,modulus=10,level=6,tiles_per_request=4,endpoint='https://gea.esac.esa.int/tap-server/tap/sync'):
    if not (0<=ra<360 and -90<=dec<=90 and 0<radius<=10 and g_min<g_max and 1<=level<=12 and tiles_per_request>=1 and modulus>=1):
        raise ValueError('Invalid sky, magnitude, tiling or sampling parameters')
    out=Path(out); meta_path=out.with_suffix('.json')
    if out.exists() or meta_path.exists(): raise ValueError('Final output already exists')
    parts=out.parent/(out.name+'.parts'); parts.mkdir(parents=True,exist_ok=True)
    T=load_study(DEFAULT_STUDY,True)
    pixels=np.sort(T.skymath.cone_tiles(level,ra,dec,radius))
    frames=[]; queries=[]; hashes={}
    for i,start in enumerate(range(0,len(pixels),tiles_per_request)):
        group=pixels[start:start+tiles_per_request]
        q=tile_query(group,level,g_min,g_max,modulus); queries.append(q)
        path=parts/f'{i:04d}.csv.gz'; metadata=path.with_suffix('.json')
        if path.exists() and metadata.exists():
            saved=json.loads(metadata.read_text(encoding='utf-8'))
            if saved['query_adql']!=q or saved['sha256']!=sha256_file(path): raise ValueError('Cached tile query/hash mismatch')
        else:
            # A partial file is diagnostic only. It must never become a cached sample.
            partial=path.with_suffix(path.suffix+'.download.vot')
            if partial.exists(): partial.unlink()
            for attempt in range(3):
                try:
                    download_parent(q,path,100000,endpoint)
                    break
                except (TimeoutError,OSError,EOFError,IncompleteRead) as exc:
                    if attempt==2: raise
                    print(f'Patch {i+1}: transient transfer failure ({type(exc).__name__}); retrying',flush=True)
                    time.sleep(5)
        frame=pd.read_csv(path,low_memory=False)
        # Astropy's older default prefers FIELD ID over name; AIP gives the
        # source_id FIELD the ID 'datalinkid'. Preserve the cached raw export,
        # recover its scientific name and verify every ID against its sky tile.
        if 'source_id' not in frame and 'datalinkid' in frame:
            frame=frame.rename(columns={'datalinkid':'source_id'})
        fac=2**35*4**(12-level)
        if 'source_id' not in frame or not pd.api.types.is_integer_dtype(frame.source_id) or not (frame.source_id//fac).isin(group).all():
            raise ValueError('Source IDs missing, inexact, or outside the requested HEALPix tiles')
        frames.append(frame)
        hashes[path.name]=sha256_file(path)
        print(f'Parent patch {i+1}/{int(np.ceil(len(pixels)/tiles_per_request))}: {len(frame)} sources',flush=True)
    merged=pd.concat(frames,ignore_index=True)
    if merged.source_id.duplicated().any(): raise ValueError('Duplicate source IDs across disjoint patches')
    a,d,a0,d0=map(np.radians,(merged.ra,merged.dec,ra,dec))
    sep=np.degrees(np.arccos(np.clip(np.sin(d0)*np.sin(d)+np.cos(d0)*np.cos(d)*np.cos(a-a0),-1,1)))
    merged=merged.loc[sep<=radius].copy()
    merged.to_csv(out,index=False,compression='infer')
    meta=dict(parent_is_unfiltered=True,synthetic_calibration=False,query_status='OK',
        source_description='Gaia DR3 disjoint indexed sky tiles with only G and independent random_index sampling; exact cone retained locally; no member or quality cuts',
        query_adql=tile_query(pixels,level,g_min,g_max,modulus),patch_queries=queries,endpoint=endpoint,
        center=[ra,dec],radius_deg=radius,g_min=g_min,g_max=g_max,random_sampling_modulus=modulus,
        random_sampling_fraction=1/modulus,n_sources=len(merged),sha256=sha256_file(out),patch_sha256=hashes)
    meta['column_name_policy']='Use VOTable FIELD names; legacy AIP datalinkid IDs recovered as source_id and checked against each requested HEALPix tile.'
    meta_path.write_text(json.dumps(meta,indent=2),encoding='utf-8')
    return meta


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ra','dec','radius'): p.add_argument('--'+k,type=float,required=True)
    p.add_argument('--out',type=Path,required=True); p.add_argument('--sampling-modulus',type=int,default=10)
    p.add_argument('--tiles-per-request',type=int,default=4)
    p.add_argument('--endpoint',choices=('esa','aip'),default='esa')
    a=p.parse_args()
    endpoint='https://gaia.aip.de/tap/sync' if a.endpoint=='aip' else 'https://gea.esac.esa.int/tap-server/tap/sync'
    result=fetch(a.out,a.ra,a.dec,a.radius,modulus=a.sampling_modulus,tiles_per_request=a.tiles_per_request,endpoint=endpoint)
    print(json.dumps(dict(out=str(a.out),n_sources=result['n_sources'],sha256=result['sha256'],endpoint=endpoint),indent=2))
