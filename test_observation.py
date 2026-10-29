import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd

from population_model import Population, build_population, read_isochrone
from response_model import ErrorLibrary, astrometric_covariance, expected_template, history_weights
from selection_bridge import DEFAULT_STUDY, load_study
from selection_model import Completeness, QualityMap, calibrate_quality, tcg_probability, parent_query, assert_tap_complete
from fit_rate import analyse, analyse_ensemble, ensemble_objective
from compare_models import compare


def parent(n=60):
    c=.8
    return pd.DataFrame({k:np.full(n,v) for k,v in dict(ra=10.,dec=0.,pmra=1.,pmdec=2.,parallax=.1,
        parallax_error=0.,pmra_error=0.,pmdec_error=0.,phot_g_mean_mag=18.,phot_bp_mean_mag=18.3,phot_rp_mean_mag=17.5,
        bp_rp=c,ruwe=1.,phot_bp_rp_excess_factor=1.162004+.011464*c+.049255*c*c-.005879*c**3,
        parallax_pmra_corr=0.,parallax_pmdec_corr=0.,pmra_pmdec_corr=0.,
        phot_g_mean_mag_error=0.,phot_bp_mean_mag_error=0.,phot_rp_mean_mag_error=0.).items()})


class ObservationChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.T=load_study(DEFAULT_STUDY,True)

    def test_population_mass_normalization_and_parsec_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=Path(tmp)/'iso.dat'
            f.write_text('# Mini Mass int_IMF Gmag G_BPmag G_RPmag\n0.5 0.5 0.0 4 4.3 3.5\n1.0 0.8 1.0 3 3.3 2.5\n',encoding='utf-8')
            pop=build_population(read_isochrone(f),mass_basis='isochrone_current_stellar_mass',source_description='SYNTHETIC',age_years=1e10,metallicity_dex=-1,synthetic=True)
            self.assertAlmostEqual(np.sum(pop.table.current_mass_msun*pop.table.number_per_msun),1.)
            pop.save(Path(tmp)/'pop')
            loaded=Population.load(Path(tmp)/'pop')
            np.testing.assert_allclose(loaded.table.number_per_msun,pop.table.number_per_msun)
            f2=Path(tmp)/'pop/population.csv'; f2.write_text(f2.read_text()+'\n',encoding='utf-8')
            with self.assertRaises(ValueError): Population.load(Path(tmp)/'pop')

    def test_invalid_population_and_age_mismatch(self):
        t=pd.DataFrame(dict(Mini=[.5,1.],Mass=[.5,.8],int_IMF=[0.,1.],Gmag=[4.,3.],G_BPmag=[4.3,3.3],G_RPmag=[3.5,2.5],logAge=[9.,9.]))
        with self.assertRaises(ValueError):
            build_population(t,mass_basis='initial_formed_mass',source_description='test',age_years=1e10,metallicity_dex=-1)

    def test_quality_counts_include_missing_sources(self):
        d=parent(); d.loc[:9,'bp_rp']=np.nan; d.loc[10:19,'ruwe']=2.
        axes=[[-1,1],[-1,1],[17,19],[0,1]]
        meta=dict(parent_is_unfiltered=True,source_description='synthetic',synthetic_calibration=True)
        q=calibrate_quality(d,self.T,center=(10,0),axes=axes,parent_metadata=meta,min_parent=20)
        pa,pq=q.query([0],[0],[18],[.5])
        self.assertAlmostEqual(pa[0],51/62)
        self.assertAlmostEqual(pq[0],41/52)
        sampled=QualityMap(q.arrays,q.metadata,draw_seed=1)
        a,b=sampled.query([0,0],[0,0],[18,18],[.5,.5])
        self.assertEqual(a[0],a[1]); self.assertEqual(b[0],b[1])
        with self.assertRaises(ValueError): q.query([5],[0],[18],[.5])
        with self.assertRaises(ValueError):
            calibrate_quality(d,self.T,center=(10,0),axes=axes,parent_metadata=dict(parent_is_unfiltered=False),min_parent=20)

    def test_covariance_and_rejection(self):
        d=parent(2); d['parallax_error']=.1; d['pmra_error']=.2; d['pmdec_error']=.3
        d['parallax_pmra_corr']=.5
        root=astrometric_covariance(d)
        C=root@np.swapaxes(root,1,2)
        np.testing.assert_allclose(C[:,0,1],.01,atol=1e-12)
        np.testing.assert_allclose(np.diagonal(C,axis1=1,axis2=2),[[.01,.04,.09]]*2,atol=1e-12)
        d['parallax_pmdec_corr']=.99; d['pmra_pmdec_corr']=-.99
        with self.assertRaises(ValueError): astrometric_covariance(d)

    def test_tcg_monotonicity_and_actual_fits_lookup(self):
        p=tcg_probability(np.array([10.,19.,21.,30.]),20.)
        self.assertTrue((np.diff(p)<0).all()); self.assertGreater(p[0],.999); self.assertLess(p[-1],1e-10)
        from astropy.table import Table
        from astropy_healpix import HEALPix
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'map.fits'
            Table(dict(UNIQ=np.arange(4,16,dtype=np.int64),M10=np.arange(12)*.1+19)).write(path)
            sf=Completeness(dict(backend='tcg_map',path=str(path),mode='multi',source_description='synthetic HEALPix layout'),tmp)
            lon,lat=HEALPix(nside=1,order='nested').healpix_to_lonlat(np.arange(12))
            got=sf.query(lon.deg,lat.deg,np.full(12,20.),np.zeros(12),np.zeros(12))
            np.testing.assert_allclose(got,tcg_probability(20.,np.arange(12)*.1+19),atol=1e-12)

    def test_rate_history_weights_and_coverage(self):
        particles=pd.DataFrame(dict(lookback_release_myr=[25.,75.],mass_weight_msun=[50.,50.]))
        run=dict(config=dict(duration_myr=100.,reference_mdot_msun_myr=1.))
        history=pd.DataFrame(dict(lookback_myr=[0.,50.,100.],escaped_stellar_rate_msun_myr=[1.,2.,3.]))
        w,ref,meta=history_weights(particles,run,history)
        np.testing.assert_allclose(w,[75,125]); self.assertEqual(ref,2.)
        with self.assertRaises(ValueError): history_weights(particles,run,history.iloc[:2])

    def test_expected_counts_zero_noise_no_double_quality_cut(self):
        d=parent(); axes=[[-1,1],[-1,1],[17,19],[0,1]]
        q=calibrate_quality(d,self.T,center=(10,0),axes=axes,
                            parent_metadata=dict(parent_is_unfiltered=True,synthetic_calibration=True,source_description='test'),min_parent=20)
        errors=ErrorLibrary(d,self.T,dict(parent_is_unfiltered=True,source_description='test'),{},(10,0))
        pop=Population(pd.DataFrame(dict(M_G=[3.],M_BP=[3.3],M_RP=[2.5],current_mass_msun=[.5],number_per_msun=[2.])),
                       dict(mass_basis='isochrone_current_stellar_mass',source_description='synthetic'))
        particles=pd.DataFrame(dict(ra_deg=[10.],dec_deg=[0.],distance_kpc=[10.],pmra_masyr=[1.],pmdec_masyr=[2.]))
        row=SimpleNamespace(ra=10.,dec=0.,R_sun=10.,pmra=1.,pmdec=2.)
        P=dict(sig0=5.,drift={'corrected':False},mids=np.array([17,18,19]),med=np.array([.8]*3),wid=np.array([.2]*3),masks=[])
        meta=dict(radius_deg=1.,pm0=[1.,2.],pm_box=2.,g_limit=20.5)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'p.npz'
            np.savez(path,x_edges=[-1,1],y_edges=[-1,1],g_edges=[17,19],probability=np.full((1,1,1),.9))
            sf=Completeness(dict(backend='tabulated',path=str(path),source_description='test'),tmp)
            signal,se,flow=expected_template(particles,pop,q,sf,errors,self.T,row,P,meta,np.ones((10,10),bool),np.linspace(-1,1,11),
                                             a_g=0,a_bp=0,a_rp=0,weights=np.array([10.]),draws=4)
            expected=20*.9*(61/62)**2
            self.assertAlmostEqual(signal.sum(),expected)
            self.assertAlmostEqual(flow['after_science'],expected)
            self.assertAlmostEqual(se.sum(),0.)

    def test_shared_ensemble_likelihood_and_identical_draws(self):
        s=np.array([1.,1.,0.,0.]); b=np.ones(4); n=np.array([110,110,10,10])
        direct=analyse(n,s,b)
        mixed=analyse_ensemble(n,np.stack([s,s]),b)
        self.assertAlmostEqual(mixed['amplitude'],direct['amplitude'],places=4)
        np.testing.assert_allclose(mixed['amplitude_interval_68_asymptotic'],direct['amplitude_interval_68_asymptotic'],atol=1e-4)
        draws=np.array([[4.,.1,.1],[.1,4.,.1]])
        objective,*_=ensemble_objective([4,0,0],draws,[1,1,1])
        lam=draws+.2
        likelihoods=np.exp(np.sum(np.array([4,0,0])*np.log(lam)-lam,axis=1))
        self.assertAlmostEqual(objective([1,.2]),-np.log(likelihoods.mean()))

    def test_predeclared_model_ranking_and_observation_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp)
            n=np.array([40.,20.,10.,10.]); b=np.full(4,5.)
            desc=dict(selection_sha256='abc',observation_identity='same_pixels',mass_basis='isochrone_current_stellar_mass')
            for name,s in (('good',[3,1,0,0]),('bad',[0,0,3,1])):
                np.savez(tmp/(name+'.npz'),counts=n,signal=s,background=b,reference_mdot=1.,
                         selection_description=json.dumps(desc),synthetic_calibration=True)
            result=compare([tmp/'good.npz',tmp/'bad.npz'],tmp/'rank')
            self.assertTrue(result['models'][0]['template'].endswith('good.npz'))
            self.assertGreater(result['models'][1]['delta_twice_nll'],0)
            desc['observation_identity']='different_pixels'
            np.savez(tmp/'changed.npz',counts=n,signal=[0,0,3,1],background=b,reference_mdot=1.,selection_description=json.dumps(desc))
            with self.assertRaises(ValueError): compare([tmp/'good.npz',tmp/'changed.npz'],tmp/'reject')

    def test_parent_query_and_tap_overflow_guard(self):
        query=parent_query(10.,0.,1.)
        where=query.split('WHERE')[1]
        self.assertNotIn('pmra',where); self.assertNotIn('ruwe',where)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'result.vot'
            path.write_text('<VOTABLE><RESOURCE><INFO name="QUERY_STATUS" value="OK"/></RESOURCE></VOTABLE>')
            assert_tap_complete(path)
            path.write_text('<VOTABLE><RESOURCE><INFO name="QUERY_STATUS" value="OK"/><INFO name="QUERY_STATUS" value="OVERFLOW"/></RESOURCE></VOTABLE>')
            with self.assertRaises(RuntimeError): assert_tap_complete(path)


if __name__=='__main__':
    unittest.main()
