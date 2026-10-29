import unittest

import numpy as np
import pandas as pd

from fit_rate import analyse, null_bootstrap
from model_extent import weighted_quantile
from tail_extent import measure_maps, project_track


class ExtentChecks(unittest.TestCase):
    def setUp(self):
        self.edges=np.linspace(-3,3,121)
        c=(self.edges[:-1]+self.edges[1:])/2
        self.X,self.Y=np.meshgrid(c,c,indexing='ij')
        self.bg=np.ones(self.X.shape)
        self.mask=np.ones(self.X.shape,bool)
        self.track=pd.DataFrame(dict(x_deg=[-4.,0.,4.],y_deg=[0.,0.,0.],side=[-1,0,1]))

    def measure(self,counts,mask=None):
        return measure_maps(counts,self.bg,self.mask if mask is None else mask,self.edges,self.track,
                            inner=.3,outer=2.8,bin_deg=.5,half_width=.2)

    def test_injected_extent_and_width(self):
        H=self.bg.copy()
        H[(np.abs(self.X)>=.3)&(np.abs(self.X)<1.8)&(np.abs(self.Y)<.1)]+=30
        profile,result=self.measure(H)
        for arm in result.values():
            self.assertAlmostEqual(arm['connected_extent_deg'],1.8)
            self.assertEqual(arm['limiting_reason'],'below_conditional_threshold')
        self.assertTrue((profile.loc[profile.candidate,'aperture_rms_width_deg']<.11).all())

    def test_no_excess_is_not_reported_as_a_tail(self):
        _,result=self.measure(self.bg)
        self.assertTrue(all(v['connected_extent_deg'] is None for v in result.values()))

    def test_mask_gap_stops_connected_extent(self):
        H=self.bg.copy()
        H[(self.X>=.3)&(self.X<2.3)&(np.abs(self.Y)<.1)]+=30
        mask=self.mask.copy()
        mask[(self.X>=1.3)&(self.X<1.8)]=False
        _,result=self.measure(H,mask)
        self.assertAlmostEqual(result['forward']['connected_extent_deg'],1.3)
        self.assertEqual(result['forward']['limiting_reason'],'coverage_or_track_limit')

    def test_projection_on_curved_polyline(self):
        v=np.array([[-1.,0.],[0.,0.],[1.,1.]])
        s=np.array([-1.,0.,np.sqrt(2)])
        along,cross,distance=project_track(np.array([[.5,.5],[-.5,.1]]),v,s,.02)
        np.testing.assert_allclose(along,[np.sqrt(.5),-.5],atol=1e-8)
        np.testing.assert_allclose(distance,[0,.1],atol=1e-8)

    def test_weighted_radial_quantile(self):
        self.assertEqual(weighted_quantile([100,1,2],[0,9,1],.9),1)
        self.assertEqual(weighted_quantile([1,2],[9,1],1),2)
        with self.assertRaises(ValueError): weighted_quantile([1],[0],.95)

    def test_profile_likelihood_intervals_and_zero_boundary(self):
        s=np.array([1.,1.,0.,0.]); b=np.ones(4)
        r=analyse([110,110,10,10],s,b)
        self.assertAlmostEqual(r['amplitude'],100,places=3)
        self.assertLess(r['amplitude_interval_68_asymptotic'][0],100)
        self.assertGreater(r['amplitude_interval_68_asymptotic'][1],100)
        self.assertGreater(r['amplitude_upper_95_one_sided_asymptotic'],r['amplitude_interval_68_asymptotic'][1])
        zero=analyse([0,0,0,0],s,b)
        self.assertEqual(zero['amplitude_interval_68_asymptotic'][0],0)
        self.assertAlmostEqual(zero['amplitude_upper_95_one_sided_asymptotic'],1.352771727/2,places=5)
        null=null_bootstrap([10,10,10,10],s,b,trials=20)
        self.assertGreater(null['null_p_conditional'],.5)


if __name__=='__main__':
    unittest.main()
