import unittest
import numpy as np
import pandas as pd
from selection_model import QualityMap
from quality_repair import pool_sparse_spatial_cells,SpatialPooledQualityMap,validate_pooling


def fixture(blue=8):
    n=np.full((2,2,1,2),100,dtype=np.int64);n[:,:,:,0]=blue
    arrays=dict(x_edges=np.array([-1.,0.,1.]),y_edges=np.array([-1.,0.,1.]),
        g_edges=np.array([17.,19.]),c_edges=np.array([-1.,.5,2.]),
        n_parent=np.full((2,2,1),200,dtype=np.int64),k_available=np.full((2,2,1),150,dtype=np.int64),
        n_available=n,k_quality=n//2)
    return QualityMap(arrays,dict(min_parent=20,selection_sha256='test',synthetic_calibration=True))


class QualityRepairTests(unittest.TestCase):
    def test_partition_keeps_supported_red_cells_and_pools_blue(self):
        q=fixture();p=pool_sparse_spatial_cells(q)
        self.assertEqual(len(p.group_n),5)
        self.assertEqual(p.group_n[p.group_ids[0,0,0,0]],32)
        self.assertEqual(len(np.unique(p.group_ids[:,:,:,1])),4)
        np.testing.assert_array_equal(p.arrays['n_parent'],q.arrays['n_parent'])

    def test_shared_draw_for_cells_in_same_group(self):
        p=pool_sparse_spatial_cells(fixture());draw=SpatialPooledQualityMap(p.arrays,p.metadata,17)
        _,prob=draw.query([-.5,.5],[-.5,.5],[18.,18.],[0.,0.])
        self.assertEqual(prob[0],prob[1])
        other=SpatialPooledQualityMap(p.arrays,p.metadata,17)
        np.testing.assert_array_equal(draw.group_probability,other.group_probability)

    def test_no_lowering_count_threshold_or_extending_support(self):
        with self.assertRaisesRegex(ValueError,'obtain more'):pool_sparse_spatial_cells(fixture(blue=1))
        p=pool_sparse_spatial_cells(fixture())
        with self.assertRaisesRegex(ValueError,'outside calibrated'):p.query([2.],[0.],[18.],[0.])

    def test_group_definitions_ignore_quality_successes(self):
        a=fixture();b=fixture();b.arrays['k_quality'][:]=0
        x=pool_sparse_spatial_cells(a);y=pool_sparse_spatial_cells(b)
        np.testing.assert_array_equal(x.group_ids,y.group_ids)
        self.assertNotEqual(x.group_probability[0],y.group_probability[0])

    def test_overlapping_or_corrupt_group_counts_are_rejected(self):
        p=pool_sparse_spatial_cells(fixture());a={k:v.copy() for k,v in p.arrays.items()};a['quality_group_n'][0]+=1
        with self.assertRaisesRegex(ValueError,'disagree'):SpatialPooledQualityMap(a,p.metadata)

    def test_holdout_successes_are_removed_before_fitting(self):
        rows=[dict(stage='quality',cell=f'{x},{y},0,0',holdout_n=2,holdout_k=0)
            for x in range(2) for y in range(2)]
        _,validation=validate_pooling(fixture(),pd.DataFrame(rows))
        # Original group: n=32,k=16. Remove held-out n=8,k=0.
        self.assertAlmostEqual(validation['cells'][0]['probability'],17/26)
        self.assertEqual(validation['cells'][0]['training_n'],24)

    def test_full_counts_cannot_hide_insufficient_training_support(self):
        rows=[dict(stage='quality',cell=f'{x},{y},0,0',holdout_n=4,holdout_k=2)
            for x in range(2) for y in range(2)]
        with self.assertRaisesRegex(ValueError,'training support'):validate_pooling(fixture(),pd.DataFrame(rows))


if __name__=='__main__':unittest.main()
