import unittest
import numpy as np
from joint_grid import aggregate_pixels
from production_study import assess,regular_axis

class ProductionTests(unittest.TestCase):
    def test_dense_axes_and_invalid_spacing(self):
        self.assertEqual(len(regular_axis(2,32,1)),31)
        self.assertEqual(len(regular_axis(100,2000,100)),20)
        self.assertEqual(regular_axis(2,3,.5),[2.,2.5,3.])
        for args in [(2,3,.3),(2,3,0),(3,2,1),(0,2,1)]:
            with self.assertRaises(ValueError):regular_axis(*args)
    def test_disjoint_count_and_state_sums(self):
        a=np.array([[1,2,3,4],[2,4,6,8.]])
        result=aggregate_pixels(a,np.array([9,9,20,20]))
        np.testing.assert_array_equal(result,[[3,7],[6,14]])
        np.testing.assert_array_equal(result.sum(axis=1),a.sum(axis=1))
    def test_invalid_mapping(self):
        with self.assertRaises(ValueError):aggregate_pixels([1,2],[-1,0])
        with self.assertRaises(ValueError):aggregate_pixels([1,2],[0])
    def test_edge_is_not_certification(self):
        rows=[dict(rate_msun_myr=r,duration_myr=t,negative_log_likelihood=r+t)
              for r in (2,4,6) for t in (250,500)]
        report=assess(dict(candidates=rows))
        self.assertTrue(all(report['best_on_boundary'].values()))
        self.assertFalse(report['physical_rate_measurement_validated'])

if __name__=='__main__':unittest.main()
