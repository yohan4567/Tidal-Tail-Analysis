"""Build a conditional science-input bundle from an official isochrone and uncut Gaia parent.

No simulated quality flags, photometry, or error distributions are generated here.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from population_model import build_population, read_isochrone
from provenance import sha256_file
from selection_bridge import DEFAULT_STUDY, load_study
from selection_model import calibrate_quality, normalize_gaia, science_available, quality_without_magnitude


def validate_quality(training, validation, quality, T, center):
    """Held-out beta-binomial diagnostics, conditional on the prespecified grid."""
    d=normalize_gaia(validation).reset_index(drop=True)
    x,y=T.gnomonic(d.ra.to_numpy(),d.dec.to_numpy(),*center)
    available=science_available(d)
    passing=np.zeros(len(d),bool)
    good=quality_without_magnitude(T,d.loc[available])
    passing[good.index.to_numpy()]=True
    xyz=np.column_stack([x,y,d.phot_g_mean_mag])
    xyzc=np.column_stack([xyz,d.g_rp])
    rows=[]; summaries={}
    for name,features,eligible,success,axes,nkey,kkey in (
        ('available',xyz,np.ones(len(d),bool),available,quality.axes[:3],'n_parent','k_available'),
        ('quality',xyzc,available,passing,quality.axes,'n_available','k_quality')):
        n=np.histogramdd(features[eligible],bins=axes)[0]
        k=np.histogramdd(features[success],bins=axes)[0]
        nt=quality.arrays[nkey]; kt=quality.arrays[kkey]
        p=(kt+1)/(nt+2)
        # Full predictive variance includes finite training-sample uncertainty.
        variance=n*p*(1-p)*(n+nt+2)/(nt+3)
        supported=(nt>=quality.minimum)&(n>0)
        z=np.divide(k-n*p,np.sqrt(variance),out=np.zeros_like(n),where=variance>0)
        summaries[name]=dict(n_holdout=int(eligible.sum()),n_in_grid=int(n.sum()),
            n_supported=int(n[supported].sum()),supported_fraction=float(n[supported].sum()/max(1,eligible.sum())),
            n_supported_cells=int(supported.sum()),
            brier_score_supported=float(np.sum((k*(1-p)**2+(n-k)*p**2)[supported])/max(1,n[supported].sum())),
            observed_success_fraction=float(k[supported].sum()/max(1,n[supported].sum())),
            predicted_success_fraction=float((n*p)[supported].sum()/max(1,n[supported].sum())),
            max_abs_predictive_z=float(np.max(np.abs(z[supported]))) if supported.any() else None)
        for index in np.argwhere(n>0):
            ix=tuple(index)
            rows.append(dict(stage=name,cell=','.join(map(str,index)),train_n=int(nt[ix]),train_k=int(kt[ix]),
                holdout_n=int(n[ix]),holdout_k=int(k[ix]),prediction=float(p[ix]),
                predictive_z=float(z[ix]),supported=bool(supported[ix])))
    return summaries,pd.DataFrame(rows)


def spatial_parent(parent,T,report):
    """Exclude crowded cluster cores outside the science fit, using sky alone."""
    x,y=T.gnomonic(parent.ra.to_numpy(),parent.dec.to_numpy(),report['center_ra_deg'],report['center_dec_deg'])
    keep=np.hypot(x,y)>1.5*report['tidal_radius_deg']
    for _,mx,my,mr in report['neighbour_masks']:
        keep &= np.hypot(x-mx,y-my)>mr
    return parent.loc[keep].copy()


def prepare(isochrone,parent,out,*,selection,study=DEFAULT_STUDY,center=(229.02,-.121),radius=5.05,ebv=.03):
    isochrone,parent,out=map(Path,(isochrone,parent,out))
    imeta=json.loads(isochrone.with_suffix('.json').read_text(encoding='utf-8'))
    pmeta=json.loads(parent.with_suffix('.json').read_text(encoding='utf-8'))
    if imeta.get('synthetic_calibration') or pmeta.get('synthetic_calibration'):
        raise ValueError('Science preparation does not accept synthetic test fixtures')
    if imeta.get('sha256')!=sha256_file(isochrone) or pmeta.get('sha256')!=sha256_file(parent):
        raise ValueError('Input checksum mismatch')
    if imeta.get('av')!=0 or pmeta.get('query_status')!='OK':
        raise ValueError('Need an unextinguished isochrone and a complete parent download')
    if not np.isfinite(ebv) or ebv<0: raise ValueError('Nonnegative E(B-V) required')
    table=read_isochrone(isochrone)
    # CMD's published Gaia EDR3 G2V coefficients, Cardelli+O'Donnell Rv=3.1.
    # This is an explicit low-extinction baseline, not a 3-D dust-map measurement.
    extinction=dict(a_g=.83627*3.1*ebv,a_bp=1.08337*3.1*ebv,a_rp=.63439*3.1*ebv)
    T=load_study(study,True)
    d=pd.read_csv(parent,low_memory=False)
    selection=Path(selection)
    report=json.loads((selection/'selection.json').read_text(encoding='utf-8'))
    if not np.allclose(center,[report['center_ra_deg'],report['center_dec_deg']],rtol=0,atol=1e-8):
        raise ValueError('Parent calibration and frozen selection centres differ')
    if report['selection_file_sha256']!=sha256_file(Path(T.HERE)/'step3_tail_search.py'):
        raise ValueError('Selection code changed; regenerate the frozen selection')
    if pmeta.get('radius_deg',0)<report['field_radius_deg']:
        raise ValueError('Parent download does not cover the science field')
    out.mkdir(parents=True,exist_ok=False)
    n_before=len(d); d=spatial_parent(d,T,report)
    donor_path=out/'parent_spatial_calibration.csv.gz'
    d.to_csv(donor_path,index=False,compression='gzip')
    pmeta=dict(pmeta,n_sources=len(d),sha256=sha256_file(donor_path),original_parent_sha256=sha256_file(parent),
        original_parent_n=n_before,client_spatial_exclusions=dict(core_radius_deg=1.5*report['tidal_radius_deg'],neighbour_masks=report['neighbour_masks']),
        source_description=pmeta['source_description']+'; sky-only core/neighbour exclusions outside science fit; no quality/member-variable cuts')
    donor_path.with_suffix('.json').write_text(json.dumps(pmeta,indent=2),encoding='utf-8')
    axes=[np.linspace(-radius,radius,3),np.linspace(-radius,radius,3),
          np.array([12.,18.,19.,20.,20.5,21.]),np.array([-5.,.5,5.])]
    # Partition before quality estimation. Never split stars by their observed flags.
    validation=np.random.default_rng(42).random(len(d))<.2
    train=d.loc[~validation].copy(); test=d.loc[validation].copy()
    qtrain=calibrate_quality(train,T,center=center,axes=axes,parent_metadata=pmeta,min_parent=20)
    stats,cells=validate_quality(train,test,qtrain,T,center)
    cells.to_csv(out/'quality_holdout_cells.csv',index=False)
    (out/'quality_validation.json').write_text(json.dumps(dict(seed=42,training_size=len(train),validation_size=len(test),
        metrics=stats,interpretation='Held-out field-star calibration diagnostic, not independent validation for the escaping stellar population.'),indent=2),encoding='utf-8')
    # Refit with all parent sources only after the held-out diagnostics have been recorded.
    quality=calibrate_quality(d,T,center=center,axes=axes,parent_metadata=pmeta,min_parent=20)
    quality.metadata.update(parent_sha256=sha256_file(donor_path),validation=stats)
    quality.save(out/'quality')
    variants=[]
    for tilt in (-1.,0.,1.):
        name=f'population_tilt_{tilt:+.0f}'.replace('+','plus').replace('-','minus')
        pop=build_population(table,mass_basis='isochrone_current_stellar_mass',weight_mode='int_imf',
            source_description=imeta['source_description'],age_years=imeta['age_years'],
            metallicity_dex=imeta['metallicity_mh'],escape_mass_tilt=tilt)
        pop.metadata.update(input_sha256=sha256_file(isochrone),official_model=imeta,
            rate_scope='Escaping represented living stars only; excludes unmodelled remnants/binaries and masses outside the isochrone.',
            mass_function_status='Assumed sensitivity family, not measured escape mass function')
        pop.save(out/name)
        # Distinguish a zero-correction sensitivity baseline from a calibrated zero-point model.
        for zp in (-.03,-.017,0.):
            key=f'tilt_{tilt:+.0f}_zp_{zp:+.3f}'.replace('+','plus').replace('-','minus')
            spec=dict(schema_version=2,population=name,quality='quality',mass_basis=pop.metadata['mass_basis'],
                completeness=dict(backend='tcg_map',path=str((Path(__file__).resolve().parent/'outputs/cache/gaia/allsky_uniq_10.fits').resolve()),mode='multi',
                    source_description='Official GaiaUnlimited TCG multi-resolution map'),
                errors=dict(path=str(donor_path.resolve()),metadata=str(donor_path.with_suffix('.json').resolve()),covariance_policy='required',
                    feature_scales=[.5,.3,1.5,1.5],neighbours=20,max_feature_distance=2.,min_donors=20),
                extinction=extinction,extinction_description='Harris2010 E(B-V) baseline, constant across field; CMD Gaia EDR3 G2V coefficients, Cardelli+ODonnell Rv=3.1; colour/spatial dependence unmodelled',
                observation_draws=32,calibration_draws=8,seed=42,parallax_zero_point_mas=zp,
                zero_point_description='Explicit scalar nuisance sensitivity scenario; not a per-source Lindegren zero-point correction',
                cache_dir='../cache/response',source_description='Official PARSEC population and empirical unfiltered Gaia DR3 quality/error donors; conditional on stated mass-function, extinction and zero-point scenarios',
                synthetic_calibration=False,scientific_status='conditional_real_input_model_not_a_validated_total_mass_loss_measurement',
                scenario=dict(escape_mass_tilt=tilt,parallax_zero_point_mas=zp))
            path=out/(key+'.json'); path.write_text(json.dumps(spec,indent=2),encoding='utf-8'); variants.append(path.name)
    manifest=dict(population_source=imeta,parent_source=pmeta,center=center,ebv=ebv,variants=variants,
        baseline='tilt_plus0_zp_minus0.017.json',
        scientific_status='Real inputs; conditional scenario inference only. Review holdout diagnostics, CMD match and sampling convergence before physical interpretation.',
        remaining_systematics=['Escape mass function is not measured; the +/-1 tilt family is illustrative, not a confidence interval.',
            'Scalar parallax offsets do not implement the full per-source Gaia zero-point model.',
            'Constant low foreground extinction ignores differential reddening and spectral dependence.',
            'Field-star quality transfer, unresolved binaries, remnants, background and dynamical assumptions require separate checks.'])
    (out/'calibration_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    return manifest


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('isochrone','parent','selection','out'): p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--study',type=Path,default=DEFAULT_STUDY)
    p.add_argument('--ra',type=float,default=229.02); p.add_argument('--dec',type=float,default=-.121)
    p.add_argument('--radius',type=float,default=5.05); p.add_argument('--ebv',type=float,default=.03)
    a=p.parse_args()
    result=prepare(a.isochrone,a.parent,a.out,selection=a.selection,study=a.study,center=(a.ra,a.dec),radius=a.radius,ebv=a.ebv)
    print(json.dumps(dict(out=str(a.out),n_variants=len(result['variants']),baseline=result['baseline'],scientific_status=result['scientific_status']),indent=2))
