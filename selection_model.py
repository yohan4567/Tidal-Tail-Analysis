"""Spatial Gaia completeness and support-aware, conditional availability/quality maps."""
import argparse
import hashlib
from importlib.metadata import version
import json
import os
import re
import time
from pathlib import Path
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from scipy.special import expit

from provenance import sha256_file
from network import open_url
from selection_bridge import REQUIRED, DEFAULT_STUDY, load_study

HERE=Path(__file__).resolve().parent
TCG_URLS={
    'hpx7':'https://zenodo.org/records/21444145/files/allsky_M10_hpx7.hdf5',
    'multi':'https://zenodo.org/records/21444145/files/allsky_uniq_10.fits'}
PARENT_COLUMNS=('source_id,ra,dec,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,'
    'phot_g_mean_flux_over_error,phot_bp_mean_flux_over_error,phot_rp_mean_flux_over_error,'
    'phot_bp_rp_excess_factor,ruwe,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,'
    'parallax_pmra_corr,parallax_pmdec_corr,pmra_pmdec_corr,astrometric_params_solved,astrometric_matched_transits')


def normalize_gaia(df):
    d=df.copy()
    if 'bp_rp' not in d and all(k in d for k in ('phot_bp_mean_mag','phot_rp_mean_mag')):
        d['bp_rp']=d.phot_bp_mean_mag-d.phot_rp_mean_mag
    if 'g_rp' not in d and 'phot_rp_mean_mag' in d:
        d['g_rp']=d.phot_g_mean_mag-d.phot_rp_mean_mag
    for band in ('g','bp','rp'):
        target=f'phot_{band}_mean_mag_error'
        if target in d:
            continue
        snr=f'phot_{band}_mean_flux_over_error'
        flux,err=f'phot_{band}_mean_flux',f'phot_{band}_mean_flux_error'
        if snr in d:
            d[target]=(2.5/np.log(10))/d[snr].where(d[snr]>0)
        elif flux in d and err in d:
            d[target]=(2.5/np.log(10))*d[err]/d[flux].where(d[flux]>0)
    return d


def science_available(df):
    cols=[*REQUIRED,'g_rp']
    if any(k not in df for k in cols):
        raise ValueError('Quality parent requires original selection columns and G-RP (phot_rp_mean_mag)')
    valid=np.isfinite(df[cols].to_numpy(float)).all(axis=1)
    for key in ('pmra_error','pmdec_error','parallax_error'):
        valid &= df[key].to_numpy()>=0
    return valid


def quality_without_magnitude(T,df):
    """Use the original quality formula; scientific magnitude cut is a separate stage."""
    old=T.G_MAX
    try:
        T.G_MAX=np.inf
        return T.quality(df)
    finally:
        T.G_MAX=old


def grid_indices(values,edges):
    values=np.asarray(values,float)
    edges=np.asarray(edges,float)
    if edges.ndim!=1 or len(edges)<2 or not np.isfinite(edges).all() or not (np.diff(edges)>0).all():
        raise ValueError('Grid edges must be finite and strictly increasing')
    index=np.searchsorted(edges,values,side='right')-1
    good=np.isfinite(values)&(index>=0)&(index<len(edges)-1)
    return np.clip(index,0,len(edges)-2),good


def lookup_indices(x,y,g,axes,color=None):
    vals=[x,y,g]+([] if color is None else [color])
    results=[grid_indices(v,e) for v,e in zip(vals,axes)]
    return tuple(r[0] for r in results),np.logical_and.reduce([r[1] for r in results])


def tcg_probability(g,m10):
    """Stable evaluation of the published TCG sigmoid using the upstream fit constants.

    Formula/model: Cantat-Gaudin et al. 2023, A&A 669 A55.
    Constants checked against GaiaUnlimited surveyTCG.py (2026-10-27).
    This standalone evaluation does not require healpy or reproduce Gaia processing.
    """
    g,m=np.broadcast_arrays(np.asarray(g,float),np.asarray(m10,float))
    pivot=20.519369625540833
    params=((.9848761394197864,.6473155510230146,.6929084598209412),
            (-.003935382139847386,.2230529402297744,-.09331877468160235),
            (.006144107896473064,.03681705933744438,.35140564525722895))
    centre,slope,shape=[a*np.minimum(m,pivot)+c*np.maximum(m-pivot,0)+b for a,b,c in params]
    if not np.isfinite(g).all() or not np.isfinite(m).all() or (m<=0).any() or (slope<=0).any() or (shape<=0).any():
        raise ValueError('Invalid/unsupported M10 or magnitude')
    # 1 - logistic(2*(G-centre)/slope)**shape, evaluated stably.
    from scipy.special import log_expit
    return -np.expm1(shape*log_expit(2*(g-centre)/slope))


class Completeness:
    def __init__(self,spec,base=HERE):
        self.spec=dict(spec); self.base=Path(base)
        self.metadata=dict(spec)
        self.backend=spec['backend']
        self._model=None
        if self.backend in ('tabulated','tcg_map'):
            self.path=(self.base/spec['path']).resolve()
            self.metadata['input_sha256']=sha256_file(self.path)
        if self.backend=='tabulated':
            with np.load(self.path,allow_pickle=False) as d:
                self.data={k:d[k] for k in d.files}
            self.axes=[self.data[k] for k in ('x_edges','y_edges','g_edges')]
            self.p=self.data['probability']
            if self.p.shape!=tuple(len(e)-1 for e in self.axes) or not np.isfinite(self.p).all() or ((self.p<0)|(self.p>1)).any():
                raise ValueError('Invalid tabulated completeness grid')
        elif self.backend=='tcg_map':
            from astropy.table import Table
            mode=spec.get('mode','multi')
            if mode=='multi':
                table=Table.read(self.path)
                self.uniq=np.asarray(table['UNIQ'],np.int64)
                self.m10=np.asarray(table['M10'],float)
                # Decode NUNIQ without a healpy dependency.
                level=(np.floor(np.log2(self.uniq/4))/2).astype(int)
                pixel=self.uniq-4*(4**level)
                self.starts=pixel*(4**(29-level))
                self.ends=self.starts+4**(29-level)
                order=np.argsort(self.starts)
                self.starts,self.ends,self.m10=self.starts[order],self.ends[order],self.m10[order]
            elif mode=='hpx7':
                import h5py
                with h5py.File(self.path) as f:
                    data=f['data'][()]
                self.hpx_order=np.asarray(data[:,0],int)
                self.hpx_index=np.asarray(data[:,1],np.int64)
                self.m10=np.asarray(data[:,2],float)
                if np.unique(self.hpx_order).size!=1:
                    raise ValueError('hpx7 map must contain a single HEALPix order')
                order=np.argsort(self.hpx_index)
                self.hpx_index,self.m10=self.hpx_index[order],self.m10[order]
            else:
                raise ValueError('tcg_map mode must be multi or hpx7')
            self.metadata.update(model='TCG 2023',implementation='standalone published-model map reader',
                                 model_source='https://arxiv.org/abs/2208.09335')
        elif self.backend=='gaiaunlimited':
            cache=(self.base/spec.get('cache_dir','outputs/cache/gaiaunlimited')).resolve()
            os.environ['GAIAUNLIMITED_DATADIR']=str(cache)
            try:
                from gaiaunlimited.selectionfunctions import DR3SelectionFunctionTCG
            except ImportError as exc:
                raise RuntimeError('GaiaUnlimited native dependencies unavailable. Use the published tcg_map reader or install requirements-observation.txt.') from exc
            mode=spec.get('mode','multi')
            if mode not in ('hpx7','multi'):
                raise ValueError('Use a precomputed hpx7/multi map here; build patch maps externally and record provenance')
            self._model=DR3SelectionFunctionTCG(mode)
            self.metadata['gaiaunlimited_version']=version('gaiaunlimited')
        else:
            raise ValueError('Unknown completeness backend')
        if not spec.get('source_description','').strip():
            raise ValueError('Completeness source_description required')

    def query(self,ra,dec,g,x,y):
        ra,dec,g,x,y=[np.asarray(a,float) for a in (ra,dec,g,x,y)]
        if not all(a.shape==g.shape for a in (ra,dec,x,y)):
            raise ValueError('Completeness query shapes differ')
        if self.backend=='tabulated':
            ix,good=lookup_indices(x,y,g,self.axes)
            if not good.all():
                raise ValueError('Completeness query outside calibrated grid')
            p=self.p[ix]
        elif self.backend=='gaiaunlimited':
            from astropy.coordinates import SkyCoord
            p=self._model.query(SkyCoord(ra,dec,unit='deg'),g.astype(float))
        else:
            from astropy_healpix import HEALPix
            import astropy.units as u
            mode=self.spec.get('mode','multi')
            level=29 if mode=='multi' else int(self.hpx_order[0])
            pixel=HEALPix(nside=2**level,order='nested').lonlat_to_healpix(ra*u.deg,dec*u.deg)
            if mode=='multi':
                i=np.searchsorted(self.starts,pixel,side='right')-1
                good=(i>=0)&(pixel<self.ends[np.maximum(i,0)])
            else:
                i=np.searchsorted(self.hpx_index,pixel)
                good=(i<len(self.hpx_index))
                i=np.minimum(i,len(self.hpx_index)-1)
                good &= self.hpx_index[i]==pixel
            if not good.all():
                raise ValueError('TCG map has no support for requested coordinates')
            p=tcg_probability(g,self.m10[i])
        p=np.asarray(p,float)
        if p.shape!=g.shape or not np.isfinite(p).all() or ((p<0)|(p>1)).any():
            raise ValueError('Completeness returned invalid probabilities')
        return p


class QualityMap:
    def __init__(self,arrays,metadata,draw_seed=None):
        self.arrays=arrays; self.metadata=metadata
        self.axes=[arrays[k] for k in ('x_edges','y_edges','g_edges','c_edges')]
        self.minimum=int(metadata['min_parent'])
        if self.minimum<1 or metadata.get('selection_sha256') is None:
            raise ValueError('Invalid quality map metadata')
        for nkey,kkey,ndim in (('n_parent','k_available',3),('n_available','k_quality',4)):
            n,k=arrays[nkey],arrays[kkey]
            if n.shape!=tuple(len(e)-1 for e in self.axes[:ndim]) or n.shape!=k.shape or (n<0).any() or (k<0).any() or (k>n).any() or not np.equal(n,np.round(n)).all() or not np.equal(k,np.round(k)).all():
                raise ValueError('Invalid quality parent/success counts')
        rng=np.random.default_rng(draw_seed) if draw_seed is not None else None
        self.probabilities={}
        for name,nkey,kkey in (('available','n_parent','k_available'),('quality','n_available','k_quality')):
            n,k=arrays[nkey],arrays[kkey]
            # Each calibration cell is drawn once, shared by all queried particles.
            self.probabilities[name]=rng.beta(k+1,n-k+1) if rng is not None else (k+1)/(n+2)

    def query(self,x,y,g,g_rp):
        ix3,good3=lookup_indices(x,y,g,self.axes[:3])
        ix4,good4=lookup_indices(x,y,g,self.axes,color=g_rp)
        if not (good3&good4).all():
            raise ValueError('Quality query outside calibrated spatial/magnitude/colour grid')
        n3=self.arrays['n_parent'][ix3]; n4=self.arrays['n_available'][ix4]
        if ((n3<self.minimum)|(n4<self.minimum)).any():
            raise ValueError('Quality calibration unsupported/low-count cells; add parent data or coarsen bins')
        pa,pq=self.probabilities['available'][ix3],self.probabilities['quality'][ix4]
        return pa,pq

    def save(self,path):
        path=Path(path); path.mkdir(parents=True,exist_ok=False)
        np.savez_compressed(path/'quality.npz',**self.arrays)
        meta=dict(self.metadata,table_sha256=sha256_file(path/'quality.npz'))
        (path/'quality.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')

    @classmethod
    def load(cls,path,draw_seed=None):
        path=Path(path)
        meta=json.loads((path/'quality.json').read_text(encoding='utf-8'))
        if meta['table_sha256']!=sha256_file(path/'quality.npz'):
            raise ValueError('Quality map hash mismatch')
        with np.load(path/'quality.npz',allow_pickle=False) as d:
            arrays={k:d[k] for k in d.files}
        return cls(arrays,meta,draw_seed)


def calibrate_quality(parent,T,*,center,axes,parent_metadata,min_parent=20):
    if not parent_metadata.get('parent_is_unfiltered',False):
        raise ValueError('Parent must explicitly be an unfiltered Gaia parent, not a PM/CMD-selected field')
    if not parent_metadata.get('source_description','').strip():
        raise ValueError('Parent source/query description required')
    if not parent_metadata.get('synthetic_calibration',False) and not parent_metadata.get('query_adql','').strip():
        raise ValueError('Real parent calibration requires the documented unfiltered ADQL query')
    if not parent_metadata.get('synthetic_calibration',False):
        query=parent_metadata['query_adql']
        where=re.split(r'\bWHERE\b',query,flags=re.IGNORECASE)
        if len(where)!=2 or re.search(r'\b(pmra|pmdec|parallax|ruwe|bp_rp|phot_bp_rp_excess_factor|astrometric_params_solved)\b',where[1],re.IGNORECASE):
            raise ValueError('Quality parent ADQL contains member/quality selection or lacks an inspectable WHERE clause')
        tokens=set(re.findall(r'\b[A-Za-z_]\w*\b',re.sub(r"'[^']*'",'',where[1]).upper()))
        allowed={'AND','OR','BETWEEN','CONTAINS','POINT','CIRCLE','BOX','POLYGON','MOD',
                 'RA','DEC','SOURCE_ID','RANDOM_INDEX','PHOT_G_MEAN_MAG'}
        if tokens-allowed:
            raise ValueError('Parent WHERE uses unsupported cuts; only sky/G and independent random_index sampling are accepted: '+str(sorted(tokens-allowed)))
        if parent_metadata.get('n_sources',len(parent))<len(parent):
            raise ValueError('Parent metadata row count is inconsistent')
    d=normalize_gaia(parent).reset_index(drop=True)
    x,y=T.gnomonic(d.ra.to_numpy(),d.dec.to_numpy(),*center)
    available=science_available(d)
    passing=np.zeros(len(d),bool)
    # No PM/parallax/CMD member cuts; quality combines original RUWE/excess conditions.
    good=quality_without_magnitude(T,d.loc[available])
    passing[np.flatnonzero(available)[np.isin(d.loc[available].index,good.index)]]=True
    xyz=np.column_stack([x,y,d.phot_g_mean_mag.to_numpy()])
    xyzc=np.column_stack([xyz,d.g_rp.to_numpy()])
    arrays=dict(x_edges=np.asarray(axes[0]),y_edges=np.asarray(axes[1]),g_edges=np.asarray(axes[2]),c_edges=np.asarray(axes[3]),
                n_parent=np.histogramdd(xyz,bins=axes[:3])[0],
                k_available=np.histogramdd(xyz[available],bins=axes[:3])[0],
                n_available=np.histogramdd(xyzc[available],bins=axes)[0],
                k_quality=np.histogramdd(xyzc[passing],bins=axes)[0])
    for key in ('n_parent','k_available','n_available','k_quality'):
        arrays[key]=arrays[key].astype(np.int64)
    meta=dict(schema_version=2,center_ra_deg=float(center[0]),center_dec_deg=float(center[1]),min_parent=min_parent,
              selection_sha256=sha256_file(Path(T.HERE)/'step3_tail_search.py'),parent_metadata=parent_metadata,
              synthetic_calibration=bool(parent_metadata.get('synthetic_calibration',False)),
              conditions='P(science columns available | Gaia,x,y,G) * P(original quality | available,x,y,G,G-RP)',
              assumptions=['Availability uses sky/G only because colours of missing sources are unobserved; population dependence must be tested.',
                           'Quality cells include original joint RUWE/flux-excess cuts, excluding the explicit scientific G limit and member PM/parallax/CMD cuts.',
                           'Transfer from field calibration to the target population is a conditional empirical assumption.'])
    return QualityMap(arrays,meta)


def fetch_map(mode,path):
    path=Path(path)
    if path.exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.part')
    with open_url(TCG_URLS[mode],timeout=60) as response,temporary.open('xb') as f:
        while chunk:=response.read(1024*1024):
            f.write(chunk)
    temporary.rename(path)
    (path.with_suffix(path.suffix+'.json')).write_text(json.dumps(dict(url=TCG_URLS[mode],sha256=sha256_file(path),mode=mode),indent=2),encoding='utf-8')


def parent_query(ra,dec,radius,g_min=0.,g_max=21.,sampling_modulus=1):
    if not np.isfinite([ra,dec,radius,g_min,g_max]).all() or not 0<=ra<360 or not -90<=dec<=90 or not 0<radius<=10 or not g_min<g_max:
        raise ValueError('Finite coordinates, 0 < radius <=10 deg and increasing magnitude limits required')
    if int(sampling_modulus)!=sampling_modulus or sampling_modulus<1:
        raise ValueError('Positive integer random sampling modulus required')
    query=(f'SELECT {PARENT_COLUMNS} FROM gaiadr3.gaia_source '
            f"WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.10f},{dec:.10f},{radius:.10f})) "
            f'AND phot_g_mean_mag >= {g_min:.8f} AND phot_g_mean_mag < {g_max:.8f}')
    if sampling_modulus>1:
        query+=f' AND MOD(random_index,{int(sampling_modulus)}) = 0'
    return query


def assert_tap_complete(path):
    found=False
    for _,element in ET.iterparse(path,events=('end',)):
        if element.tag.rsplit('}',1)[-1]=='INFO' and element.attrib.get('name','').upper()=='QUERY_STATUS':
            found=True
            value=element.attrib.get('value','').upper()
            if value!='OK':
                raise RuntimeError(f'TAP query incomplete ({value}): {element.text}. Use a smaller patch; do not calibrate truncated parents.')
        element.clear()
    if not found:
        raise RuntimeError('No TAP QUERY_STATUS found; cannot verify complete parent response')


def download_parent(query,out,max_records=250000,endpoint='https://gea.esac.esa.int/tap-server/tap/sync'):
    """Read-only bounded TAP request; reject VOTable overflow before exporting a parent."""
    from astropy.table import Table
    out=Path(out)
    metadata_path=out.with_suffix('.json')
    temporary=out.with_suffix(out.suffix+'.download.vot')
    if out.exists() or metadata_path.exists() or temporary.exists() or max_records<1:
        raise ValueError('Choose fresh parent output paths and a positive row limit')
    out.parent.mkdir(parents=True,exist_ok=True)
    if endpoint not in ('https://gea.esac.esa.int/tap-server/tap/sync','https://gaia.aip.de/tap/sync'):
        raise ValueError('Use a supported Gaia DR3 TAP endpoint')
    body=urllib.parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='votable',MAXREC=max_records,QUERY=query)).encode()
    request=urllib.request.Request(endpoint,data=body,headers={'Content-Type':'application/x-www-form-urlencoded'})
    started=time.monotonic()
    try:
        with open_url(request,timeout=60) as response,temporary.open('xb') as f:
            read=getattr(response,'read1',response.read)
            while chunk:=read(64*1024):
                f.write(chunk)
                if time.monotonic()-started>300:
                    raise TimeoutError('Parent transfer exceeded 300 seconds; use resumable smaller sky tiles')
        assert_tap_complete(temporary)
        frame=Table.read(temporary,format='votable',use_names_over_ids=True).to_pandas()
    except Exception:
        if temporary.exists(): temporary.unlink()
        raise
    frame.columns=[str(k).lower() for k in frame.columns]
    frame=normalize_gaia(frame)
    frame.to_csv(out,index=False,compression='infer')
    meta=dict(parent_is_unfiltered=True,source_description='Gaia DR3 parent with only recorded sky/G limits; no quality or member cuts',
              query_adql=query,endpoint=endpoint,max_records=max_records,n_sources=len(frame),
              synthetic_calibration=False,sha256=sha256_file(out),query_status='OK')
    sample=re.search(r'MOD\(random_index,(\d+)\)',query,re.IGNORECASE)
    if sample:
        meta.update(random_sampling_modulus=int(sample.group(1)),random_sampling_fraction=1/int(sample.group(1)),
            source_description='Gaia DR3 sky/G parent with independent random_index subsampling; no quality or member cuts')
    metadata_path.write_text(json.dumps(meta,indent=2),encoding='utf-8')
    temporary.unlink()
    return meta


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='action',required=True)
    fetch=sub.add_parser('fetch-gaia-map'); fetch.add_argument('--mode',choices=TCG_URLS,default='multi'); fetch.add_argument('--out',type=Path,required=True)
    query=sub.add_parser('parent-query')
    for key in ('ra','dec','radius'):
        query.add_argument('--'+key,type=float,required=True)
    query.add_argument('--g-min',type=float,default=0.); query.add_argument('--g-max',type=float,default=21.)
    query.add_argument('--out',type=Path,required=True); query.add_argument('--download',action='store_true')
    query.add_argument('--max-records',type=int,default=250000)
    query.add_argument('--endpoint',choices=('esa','aip'),default='esa')
    query.add_argument('--sampling-modulus',type=int,default=1,help='Unbiased independent random_index subsampling for bounded parent downloads')
    train=sub.add_parser('calibrate-quality')
    for key in ('parent','parent-metadata','grid','out'):
        train.add_argument('--'+key,type=Path,required=True)
    train.add_argument('--study',type=Path,default=DEFAULT_STUDY)
    train.add_argument('--min-parent',type=int,default=20)
    a=p.parse_args()
    if a.action=='fetch-gaia-map':
        fetch_map(a.mode,a.out); print(a.out); return
    if a.action=='parent-query':
        text=parent_query(a.ra,a.dec,a.radius,a.g_min,a.g_max,a.sampling_modulus)
        if a.download:
            endpoint='https://gaia.aip.de/tap/sync' if a.endpoint=='aip' else 'https://gea.esac.esa.int/tap-server/tap/sync'
            print(json.dumps(download_parent(text,a.out,a.max_records,endpoint),indent=2))
        else:
            a.out.parent.mkdir(parents=True,exist_ok=True)
            with a.out.open('x',encoding='utf-8') as f: f.write(text+'\n')
            print(a.out)
        return
    meta=json.loads(a.parent_metadata.read_text(encoding='utf-8'))
    grid=json.loads(a.grid.read_text(encoding='utf-8'))
    T=load_study(a.study)
    quality=calibrate_quality(pd.read_csv(a.parent,low_memory=False),T,center=grid['center'],
                              axes=[np.asarray(grid[k],float) for k in ('x_edges','y_edges','g_edges','c_edges')],
                              parent_metadata=meta,min_parent=a.min_parent)
    quality.metadata['parent_sha256']=sha256_file(a.parent)
    quality.save(a.out)
    print(json.dumps(quality.metadata,indent=2))


if __name__=='__main__':
    main()
