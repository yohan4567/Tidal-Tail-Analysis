import unittest
from types import SimpleNamespace

import numpy as np
import pandas as pd

from observe_particles import make_mock
from selection_bridge import apply_selection, check_catalogue, load_study, DEFAULT_STUDY


def catalogue(n=30):
    c=.8
    exc=1.162004+.011464*c+.049255*c*c-.005879*c**3
    return pd.DataFrame({k:np.full(n,v) for k,v in dict(ra=10.,dec=0.,pmra=1.,pmdec=2.,
        pmra_error=0.,pmdec_error=0.,parallax=.1,parallax_error=0.,phot_g_mean_mag=18.,
        bp_rp=c,ruwe=1.,phot_bp_rp_excess_factor=exc).items()})


class SelectionChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.T=load_study(DEFAULT_STUDY,legacy=True)

    def test_observation_normalization_and_zero_error_limit(self):
        particles=pd.DataFrame(dict(ra_deg=[10.,10.],dec_deg=[0.,0.],distance_kpc=[10.,10.],
                                    pmra_masyr=[1.,1.],pmdec_masyr=[2.,2.],mass_weight_msun=[3.,3.]))
        pop=pd.DataFrame(dict(M_G=[3.],bp_rp=[.8],number_per_msun=[2.]))
        comp=pd.DataFrame(dict(G=[17.,19.],completeness=[.5,.5]))
        mock=make_mock(particles,pop,catalogue(),comp,a_g=0,reddening=0)
        np.testing.assert_allclose(mock.expected_count,[3.,3.])
        np.testing.assert_allclose(mock.parallax,.1)
        np.testing.assert_allclose(mock.pmra,1.)
        np.testing.assert_allclose(mock.phot_g_mean_mag,18.)
        self.assertEqual(len(self.T.quality(mock)),2)
        with self.assertRaises(ValueError):
            make_mock(particles,pop,catalogue(),pd.DataFrame(dict(G=[19.,20.],completeness=[1.,1.])),a_g=0,reddening=0)

    def test_frozen_cmd_download_and_neighbour_masks(self):
        d=catalogue(5)
        d.loc[1,'pmra']=10.  # outside downloaded PM window
        d.loc[2,'ra']=10.5  # neighbour mask
        d.loc[3,'bp_rp']=2. # CMD rejection (also poor excess)
        d.loc[4,'phot_g_mean_mag']=21.
        P=dict(masks=[('other',.5,0,.1)],drift={'corrected':False},sig0=5.,
               mids=np.array([17.,18.,19.]),med=np.array([.8,.8,.8]),wid=np.array([.1,.1,.1]))
        row=SimpleNamespace(ra=10.,dec=0.,pmra=1.,pmdec=2.,R_sun=10.)
        meta=dict(radius_deg=1.,pm0=[1.,2.],pm_box=2.,g_limit=20.5)
        selected,flow=apply_selection(self.T,d,row,P,meta)
        self.assertEqual(len(selected),1)
        self.assertEqual(flow['download_window'],3)
        with self.assertRaises(ValueError): check_catalogue(pd.DataFrame({'ra':[1.]}))


if __name__=='__main__':
    unittest.main()
