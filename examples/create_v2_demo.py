"""Create explicitly synthetic calibration fixtures for the hybrid observation pipeline."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from population_model import build_population
from selection_model import calibrate_quality
from selection_bridge import load_study, DEFAULT_STUDY


def create(out,study=DEFAULT_STUDY,real_gaia_map=None,variable_rate=False):
    out=Path(out).resolve(); out.mkdir(parents=True,exist_ok=False)
    rng=np.random.default_rng(42); count=8000
    ra0,dec0=229.02,-.121
    g=rng.uniform(16,21,count); c=rng.normal(.5,.025,count); bp_rp=rng.normal(1.02,.025,count)
    parent=pd.DataFrame(dict(ra=ra0+rng.uniform(-5,5,count)/np.cos(np.radians(dec0)),
        dec=dec0+rng.uniform(-5,5,count),phot_g_mean_mag=g,phot_rp_mean_mag=g-c,
        phot_bp_mean_mag=g-c+bp_rp,bp_rp=bp_rp,pmra=rng.normal(-2.754,1,count),
        pmdec=rng.normal(-2.643,1,count),parallax=rng.normal(.05,.1,count),
        pmra_error=.06+np.maximum(g-17,0)*.03,pmdec_error=.07+np.maximum(g-17,0)*.03,
        parallax_error=.05+np.maximum(g-17,0)*.02,ruwe=np.where(rng.random(count)<.8,1.,2.),
        parallax_pmra_corr=np.full(count,.1),parallax_pmdec_corr=np.full(count,-.1),pmra_pmdec_corr=np.full(count,.2),
        phot_g_mean_mag_error=np.full(count,.01),phot_bp_mean_mag_error=np.full(count,.02),phot_rp_mean_mag_error=np.full(count,.02)))
    parent['phot_bp_rp_excess_factor']=1.162004+.011464*bp_rp+.049255*bp_rp**2-.005879*bp_rp**3
    missing=rng.random(count)<.1
    parent.loc[missing,['phot_bp_mean_mag','bp_rp']]=np.nan
    parent.to_csv(out/'parent_SYNTHETIC.csv',index=False)
    meta=dict(parent_is_unfiltered=True,synthetic_calibration=True,
              source_description='SYNTHETIC software-test population with independent 90% availability and 80% quality; not actual Gaia calibration')
    (out/'parent.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    iso=pd.DataFrame(dict(initial_mass_msun=[.7,.75,.8],current_mass_msun=[.7,.75,.8],
                          M_G=[2.,2.1,2.2],M_BP=[2.52,2.62,2.72],M_RP=[1.5,1.6,1.7],int_IMF=[0.,.3,.6]))
    iso.to_csv(out/'isochrone_SYNTHETIC.csv',index=False)
    pop=build_population(iso,mass_basis='isochrone_current_stellar_mass',source_description='SYNTHETIC arbitrary isochrone for integration tests',
                         age_years=12e9,metallicity_dex=-1.4,synthetic=True)
    pop.save(out/'population')
    axes=[np.array([-6.,0.,6.]),np.array([-6.,0.,6.]),np.array([16.,17.,18.,19.,20.,21.]),np.array([0.,1.5])]
    quality=calibrate_quality(parent,load_study(study),center=(ra0,dec0),axes=axes,parent_metadata=meta,min_parent=20)
    quality.save(out/'quality')
    if real_gaia_map:
        completeness=dict(backend='tcg_map',path=str(Path(real_gaia_map).resolve()),mode='multi',
                          source_description='Official TCG reference map, with synthetic population and quality calibration elsewhere in this test')
    else:
        np.savez_compressed(out/'survey_SYNTHETIC.npz',x_edges=axes[0],y_edges=axes[1],g_edges=np.array([0.,40.]),
                            probability=np.full((2,2,1),.9))
        completeness=dict(backend='tabulated',path='survey_SYNTHETIC.npz',synthetic_calibration=True,
                          source_description='SYNTHETIC constant 0.9 completeness for software tests')
    spec=dict(schema_version=2,population='population',quality='quality',mass_basis='isochrone_current_stellar_mass',
              completeness=completeness,errors=dict(path='parent_SYNTHETIC.csv',metadata='parent.json',covariance_policy='required'),
              extinction=dict(a_g=0.,a_bp=0.,a_rp=0.),observation_draws=16,calibration_draws=4,seed=42,
              parallax_zero_point_mas=0.,zero_point_description='SYNTHETIC zero-correction baseline; not a calibrated Gaia zero point',
              cache_dir='../cache/response',
              source_description='SYNTHETIC integration test only; all astrophysical inputs arbitrary',synthetic_calibration=True)
    if variable_rate:
        pd.DataFrame(dict(lookback_myr=[0.,50.,1000.],escaped_stellar_rate_msun_myr=[1.,2.,1.])).to_csv(out/'rate_history_SYNTHETIC.csv',index=False)
        spec.update(rate_history='rate_history_SYNTHETIC.csv',rate_history_mass_basis=spec['mass_basis'],
                    rate_history_source='SYNTHETIC arbitrary time history for interface tests; not an EMACSS prediction',
                    assume_uniform_release_for_legacy=True)
    (out/'response.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
    print(out/'response.json')
    return out/'response.json'


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--study',type=Path,default=DEFAULT_STUDY)
    p.add_argument('--real-gaia-map',type=Path)
    p.add_argument('--variable-rate-demo',action='store_true')
    a=p.parse_args()
    create(a.out,a.study,a.real_gaia_map,a.variable_rate_demo)
