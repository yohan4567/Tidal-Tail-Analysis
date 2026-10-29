"""Reuse tidal_study_v2 selection on observations and Gaia-like mock catalogues."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from provenance import sha256_file

HERE = Path(__file__).resolve().parent
DEFAULT_STUDY = Path(os.environ.get('GC_TIDAL_STUDY',str(HERE.parent/'TidaleffectsOnGCs-main/TidaleffectsOnGCs-main/tidal_study_v2'))).expanduser()
REQUIRED = ('ra','dec','pmra','pmdec','pmra_error','pmdec_error','parallax',
            'parallax_error','phot_g_mean_mag','bp_rp','ruwe','phot_bp_rp_excess_factor')


def load_study(path, legacy=False):
    """Load existing functions without exposing this CLI's arguments to them."""
    os.environ.setdefault('MPLCONFIGDIR',str(HERE/'outputs/.mplcache'))
    path = Path(path).resolve()
    if not (path/'step3_tail_search.py').is_file():
        raise FileNotFoundError('External tidal_study_v2 selection code is not bundled. Set GC_TIDAL_STUDY or pass --study to its directory; see docs/SETUP_AND_DATA.md.')
    old_argv, old_path = sys.argv, sys.path.copy()
    try:
        sys.argv = ['selection_bridge'] + (['--v1'] if legacy else [])
        sys.path.insert(0,str(path))
        spec=importlib.util.spec_from_file_location('_gc_stream_tail_selection',path/'step3_tail_search.py')
        module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.argv=old_argv
        sys.path[:]=old_path
    return module


def check_catalogue(df):
    missing=sorted(set(REQUIRED)-set(df.columns))
    if missing:
        raise ValueError('Gaia-like catalogue missing columns: '+', '.join(missing)+
                         '. Raw particles.csv needs a photometric/error observation model first.')
    for k in ('pmra_error','pmdec_error','parallax_error'):
        if (df[k].dropna()<0).any():
            raise ValueError(f'Negative uncertainty: {k}')


def apply_selection(T, df, row, P, meta, missing_tiles=None, download=True):
    """Freeze real-data CMD, drift decision and masks; return cumulative cut flow."""
    check_catalogue(df)
    d=df.copy().reset_index(drop=True)
    flow={'input':len(d)}
    if download:
        if meta is None or not all(k in meta for k in ('radius_deg','pm0','pm_box','g_limit')):
            raise ValueError('Download metadata required to reproduce the observed footprint/window')
        a,z,a0,z0=map(np.radians,(d.ra,d.dec,row.ra,row.dec))
        sep=np.degrees(np.arccos(np.clip(np.sin(z0)*np.sin(z)+np.cos(z0)*np.cos(z)*np.cos(a-a0),-1,1)))
        keep=(sep<=meta['radius_deg']) & (d.pmra>=meta['pm0'][0]-meta['pm_box']) & (d.pmra<=meta['pm0'][0]+meta['pm_box'])
        keep &= (d.pmdec>=meta['pm0'][1]-meta['pm_box']) & (d.pmdec<=meta['pm0'][1]+meta['pm_box']) & (d.phot_g_mean_mag<meta['g_limit'])
        d=d[keep].copy()
    flow['download_window']=len(d)
    d=T.quality(d)
    flow['quality']=len(d)
    d['x'],d['y']=T.gnomonic(d.ra.to_numpy(),d.dec.to_numpy(),row.ra,row.dec)
    d['r']=np.hypot(d.x,d.y)
    for _,mx,my,mr in P['masks']:
        d=d[np.hypot(d.x-mx,d.y-my)>mr].copy()
    flow['neighbour_masks']=len(d)
    if P['drift']['corrected'] and len(d):
        pm,dist,_,_=T.orbit_expectation(row,P['R'],d.x.to_numpy(),d.y.to_numpy())
        dg=5*np.log10(dist/row.R_sun)
    else:
        pm,dist,dg=(row.pmra,row.pmdec),row.R_sun,None
    # Compute orbit arrays once, then apply successive masks on the same index.
    plx=np.asarray(T.plx_select(d,dist),bool)
    motion=np.asarray(T.pm_select(d,pm,dist,P['sig0']),bool)
    cmd=np.asarray(T.cmd_select(d,P['mids'],P['med'],P['wid'],dg),bool)
    flow['parallax']=int(plx.sum())
    flow['proper_motion']=int((plx & motion).sum())
    keep=plx & motion & cmd
    if missing_tiles is not None:
        tiles=T.skymath.ang2pix_nest(missing_tiles['level'],d.ra.to_numpy(),d.dec.to_numpy())
        keep &= ~np.isin(tiles,missing_tiles['tiles'])
    d=d[keep].copy()
    flow['cmd_and_missing_tiles']=len(d)
    return d,flow


def histogram(d,edges,weights=None):
    return np.histogram2d(d.x,d.y,bins=[edges,edges],weights=weights)[0]


def run(args):
    os.environ.setdefault('MPLCONFIGDIR',str(HERE/'outputs/.mplcache'))
    T=load_study(args.study,args.legacy)
    master=pd.read_csv(args.study/'catalogs/master_clusters.csv').set_index('id')
    row=master.loc[args.cluster].rename(args.cluster)
    field=T.field_file(args.cluster,'_field.csv.gz')
    if not field.exists():
        raise FileNotFoundError(f'{field}. Restore the planned field, or explicitly use --legacy for the original field.')
    metadata_path=T.field_file(args.cluster,'_field.json')
    meta=json.loads(metadata_path.read_text(encoding='utf-8'))
    P=T.prepare(args.cluster,row,master)
    mt=T.field_file(args.cluster,'_missing_tiles.json')
    missing=json.loads(mt.read_text(encoding='utf-8')) if mt.exists() else None
    raw=pd.read_csv(field,low_memory=False)
    selected,flow=apply_selection(T,raw,row,P,meta,missing,download=False)
    # Reproduce the existing sample exactly, independent of row order/source IDs.
    cols=(['source_id'] if 'source_id' in raw else [])+list(REQUIRED)
    expected=P['cand'][cols].reset_index(drop=True)
    pd.testing.assert_frame_equal(selected[cols].reset_index(drop=True),expected,check_dtype=False)
    H=histogram(selected,P['edges'])
    _,background,_=P['build'](P['cand'],P['s_pix'])
    fitmask=P['inside'] & (P['RR']>1.5*P['rt']) & (P['RR']<0.9*P['Rt'])
    if not fitmask.any():
        raise ValueError('No unmasked tail-region bins; revise field or fit region')
    out=args.out
    out.mkdir(parents=True,exist_ok=False)
    selected.to_csv(out/'observed_candidates.csv.gz',index=False,compression='gzip')
    track=T.skymath.orbit_track(row.ra,row.dec,row.R_sun,row.pmra,row.pmdec,row.rv,radius=P['R']+1.)
    tx,ty=T.gnomonic(track['ra'],track['dec'],row.ra,row.dec)
    pd.DataFrame(dict(x_deg=tx,y_deg=ty,side=track['side'],ra_deg=track['ra'],dec_deg=track['dec'],distance_kpc=track['d'])).to_csv(out/'orbit_track.csv',index=False)
    report=dict(cluster=args.cluster,legacy=args.legacy,source_field=str(field.resolve()),
                schema_version=2,sig0_kms=P['sig0'],neighbour_masks=P['masks'],missing_tiles=missing,
                cluster_selection_parameters={k:float(row[k]) for k in ('ra','dec','R_sun','pmra','pmdec','rv')},
                source_field_sha256=sha256_file(field),
                center_ra_deg=float(row.ra),center_dec_deg=float(row.dec),distance_kpc=float(row.R_sun),
                tidal_radius_deg=P['rt'],field_radius_deg=P['R'],test_radius_deg=P['Rt'],
                selection_file_sha256=hashlib.sha256((args.study/'step3_tail_search.py').read_bytes()).hexdigest(),
                skymath_sha256=hashlib.sha256((args.study/'skymath.py').read_bytes()).hexdigest(),
                source_metadata=meta,observation_cutflow=flow,exact_selection_match=True,
                drift=P['drift'],cmd_enabled=len(P['mids'])>=3,
                cmd_ridge=dict(G=P['mids'].tolist(),color=P['med'].tolist(),width=P['wid'].tolist()),
                n_fit_bins=int(fitmask.sum()),n_fit_stars=int(H[fitmask].sum()),
                fit_region='unmasked projected 1.5*r_t < R < 0.9*min(field_radius,5 deg)',
                background_caveat='Original background shape estimated from the same field; may absorb stream stars. Shape uncertainty is not propagated.',
                population_calibration='No automatic mass-to-visible-star conversion. Calibrated mock expected_count is required for mass-loss fitting.')
    signal=np.zeros_like(H)
    if args.mock is not None:
        if args.reference_mdot is None or not np.isfinite(args.reference_mdot) or args.reference_mdot<=0 or not args.calibration_description:
            raise ValueError('--mock requires positive --reference-mdot and --calibration-description')
        mock=pd.read_csv(args.mock,low_memory=False)
        sidecar=args.mock.with_suffix('.json')
        if sidecar.exists():
            provenance=json.loads(sidecar.read_text(encoding='utf-8'))
            if 'reference_mdot' in provenance and not np.isclose(provenance['reference_mdot'],args.reference_mdot,rtol=1e-12,atol=0):
                raise ValueError('Reference rate disagrees with mock provenance')
            if 'cluster' in provenance and provenance['cluster']!=args.cluster:
                raise ValueError('Mock cluster disagrees with selected observation')
            report['mock_provenance']=provenance
        if 'expected_count' not in mock or not np.isfinite(mock.expected_count).all() or (mock.expected_count<0).any():
            raise ValueError('Each mock row needs finite nonnegative expected_count at the reference rate, including population/completeness calibration')
        mock_sel,mockflow=apply_selection(T,mock,row,P,meta,missing)
        signal=histogram(mock_sel,P['edges'],mock_sel.expected_count)
        mock_sel.to_csv(out/'selected_mock.csv.gz',index=False,compression='gzip')
        description=json.dumps(dict(cluster=args.cluster,reference_mdot=args.reference_mdot,
                                    calibration=args.calibration_description,selection_sha256=report['selection_file_sha256'],
                                    source_field=str(field.resolve()),region=report['fit_region'],mock_provenance=report.get('mock_provenance')),ensure_ascii=False)
        np.savez_compressed(out/'calibrated_counts.npz',counts=H[fitmask],signal=signal[fitmask],background=background[fitmask],
                            reference_mdot=args.reference_mdot,selection_description=description,
                            synthetic_calibration=bool(report.get('mock_provenance',{}).get('synthetic_calibration',False)))
        report.update(mock_cutflow=mockflow,mock_sha256=hashlib.sha256(args.mock.read_bytes()).hexdigest(),calibration_description=args.calibration_description,
                      selected_expected_stars=float(signal[fitmask].sum()),reference_mdot=args.reference_mdot)
    np.savez_compressed(out/'maps.npz',counts=H,background=background,signal=signal,edges=P['edges'],fitmask=fitmask,inside=P['inside'])
    (out/'selection.json').write_text(json.dumps(T.clean(report),indent=2,ensure_ascii=False),encoding='utf-8')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,3,figsize=(14,4),constrained_layout=True)
    extent=[-P['R'],P['R'],-P['R'],P['R']]
    for ax,arr,title in zip(axs,[H,background,signal],['Observed raw counts','Original background estimate','Selected mock expected counts' if args.mock else 'Mock not supplied']):
        im=ax.imshow(np.where(P['inside'],arr,np.nan).T,origin='lower',extent=extent)
        ax.set(title=title,xlabel='Gnomonic east (deg)',ylabel='Gnomonic north (deg)')
        fig.colorbar(im,ax=ax,label='Stars / pixel')
    fig.savefig(out/'selection.png',dpi=150)
    plt.close(fig)
    print(json.dumps(T.clean(report),ensure_ascii=False,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,default=DEFAULT_STUDY)
    p.add_argument('--cluster',default='Pal5')
    p.add_argument('--legacy',action='store_true',help='Explicitly use original field instead of v3 field plan')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--mock',type=Path,help='Gaia-like mock catalogue with calibrated expected_count per row')
    p.add_argument('--reference-mdot',type=float,help='Total both-arm Msun/Myr represented by mock')
    p.add_argument('--calibration-description',help='Document population, detection completeness and error model')
    run(p.parse_args())


if __name__=='__main__':
    main()
