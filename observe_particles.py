"""Conditional observation model: calibrated luminosity function + empirical Gaia errors."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from provenance import sha256_file


def draw_population(particles, population, rng):
    """Number-weighted luminosity sampling; weights integrate to stars per Msun."""
    for key in ('M_G','bp_rp','number_per_msun'):
        if key not in population or not np.isfinite(population[key]).all():
            raise ValueError(f'Population requires finite {key}')
    w=population.number_per_msun.to_numpy()
    if (w<0).any() or w.sum()<=0:
        raise ValueError('Population number_per_msun must be nonnegative with positive sum')
    indices=rng.choice(len(population),size=len(particles),p=w/w.sum())
    return population.iloc[indices].reset_index(drop=True),w.sum()


def make_mock(particles,population,field,completeness,*,a_g,reddening,seed=42):
    from selection_bridge import REQUIRED, check_catalogue
    check_catalogue(field)
    required=('ra_deg','dec_deg','distance_kpc','pmra_masyr','pmdec_masyr','mass_weight_msun')
    if any(k not in particles for k in required):
        raise ValueError('Missing stream particle coordinates or mass weights')
    if not np.isfinite(particles[list(required)]).all().all() or (particles.distance_kpc<=0).any() or (particles.mass_weight_msun<0).any():
        raise ValueError('Particle values must be finite with positive distances and nonnegative mass')
    if not np.isfinite([a_g,reddening]).all() or a_g<0 or reddening<0:
        raise ValueError('Finite nonnegative extinction/reddening required')
    rng=np.random.default_rng(seed)
    pop,number=draw_population(particles,population,rng)
    g=pop.M_G.to_numpy()+5*np.log10(particles.distance_kpc.to_numpy())+10+a_g
    color=pop.bp_rp.to_numpy()+reddening
    comp=completeness.sort_values('G')
    if len(comp)<2 or not np.isfinite(comp[['G','completeness']]).all().all() or (np.diff(comp.G)<=0).any() or not comp.completeness.between(0,1).all():
        raise ValueError('Completeness requires unique G grid and probabilities in [0,1]')
    if len(g)==0 or g.min()<comp.G.min() or g.max()>comp.G.max():
        raise ValueError('Completeness table must cover every simulated apparent G; no extrapolation')
    efficiency=np.interp(g,comp.G,comp.completeness)
    # Draw joint empirical RUWE/excess/error tuples before applying quality cuts.
    donor=field.dropna(subset=list(REQUIRED)).copy()
    donor=donor[(donor.pmra_error>=0)&(donor.pmdec_error>=0)&(donor.parallax_error>=0)]
    if len(donor)<5:
        raise ValueError('Insufficient raw field error donors')
    features=np.column_stack([donor.phot_g_mean_mag.to_numpy()/0.5,donor.bp_rp.to_numpy()/0.3])
    distances,indices=cKDTree(features).query(np.column_stack([g/0.5,color/0.3]),k=min(20,len(donor)))
    # Reject extrapolation, including faint populations absent from the field.
    if np.any((efficiency>0)&(distances[:,0]>2)):
        raise ValueError('No nearby Gaia error donors for some detectable mock stars; improve field/error model')
    options=np.maximum((distances<=2).sum(axis=1),1)
    choice=(rng.random(len(g))*options).astype(int)
    donor=donor.iloc[indices[np.arange(len(g)),choice]].reset_index(drop=True)
    out=pd.DataFrame(dict(ra=particles.ra_deg.to_numpy(),dec=particles.dec_deg.to_numpy(),
                          phot_g_mean_mag=g,bp_rp=color,ruwe=donor.ruwe.to_numpy(),
                          pmra_error=donor.pmra_error.to_numpy(),pmdec_error=donor.pmdec_error.to_numpy(),
                          parallax_error=donor.parallax_error.to_numpy()))
    # Preserve excess residual around the colour-dependent quality-cut locus.
    def locus(c):
        return np.where(c<.5,1.154360+.033772*c+.032277*c*c,
                        np.where(c<4,1.162004+.011464*c+.049255*c*c-.005879*c**3,1.057572+.140537*c))
    out['phot_bp_rp_excess_factor']=locus(color)+donor.phot_bp_rp_excess_factor.to_numpy()-locus(donor.bp_rp.to_numpy())
    out['pmra']=particles.pmra_masyr.to_numpy()+rng.normal(size=len(g))*out.pmra_error
    out['pmdec']=particles.pmdec_masyr.to_numpy()+rng.normal(size=len(g))*out.pmdec_error
    out['parallax']=1/particles.distance_kpc.to_numpy()+rng.normal(size=len(g))*out.parallax_error
    out['expected_count']=particles.mass_weight_msun.to_numpy()*number*efficiency
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True,help='Simulation output directory containing run.json and particles.csv')
    p.add_argument('--population',type=Path,required=True,help='CSV: M_G, intrinsic bp_rp, number_per_msun per row')
    p.add_argument('--field',type=Path,required=True,help='Raw Gaia field BEFORE scientific selection')
    p.add_argument('--completeness',type=Path,required=True,help='CSV: apparent G, completeness (detection/catalogue inclusion, before reused quality cuts)')
    p.add_argument('--a-g',type=float,required=True)
    p.add_argument('--reddening',type=float,required=True,help='E(BP-RP)')
    p.add_argument('--calibration-description',required=True)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--synthetic-calibration',action='store_true',help='Explicitly mark software-test population/completeness inputs')
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if not a.calibration_description.strip():
        p.error('Document the luminosity function and completeness calibration')
    if a.out.exists() or a.out.with_suffix('.json').exists():
        raise FileExistsError('Choose new output paths')
    config=json.loads((a.run/'run.json').read_text(encoding='utf-8'))
    particles=pd.read_csv(a.run/'particles.csv')
    mock=make_mock(particles,pd.read_csv(a.population),pd.read_csv(a.field,low_memory=False),pd.read_csv(a.completeness),
                   a_g=a.a_g,reddening=a.reddening,seed=a.seed)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    mock.to_csv(a.out,index=False)
    report=dict(cluster=config['config']['name'],reference_mdot=config['config']['reference_mdot_msun_myr'],seed=a.seed,
                synthetic_calibration=a.synthetic_calibration,
                calibration_description=a.calibration_description,a_g=a.a_g,reddening=a.reddening,
                input_sha256={str(f.resolve()):sha256_file(f) for f in (a.population,a.field,a.completeness,a.run/'particles.csv')},
                limitations=['Independent Gaussian PM/parallax errors; no covariance or spatially correlated systematics.',
                             'Photometry deterministic per luminosity draw; no photometric scatter or spatially varying extinction.',
                             'Empirical raw-field error/quality donors matched by G and colour, not sky position.',
                             'Luminosity function must describe escaping population per total escaped mass, including unseen stars/remnants in its normalization.',
                             'One luminosity/error realization per tracer; increase tracers and repeat seeds to test convergence.',
                             'Population and completeness are external calibrations, not inferred by this tool.'])
    a.out.with_suffix('.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(f'{len(mock)} mock tracers -> {a.out}; reference rate {report["reference_mdot"]} Msun/Myr')


if __name__=='__main__':
    main()
