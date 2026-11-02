import copy
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd

from bright_calibration import combine_parent_strata, verify_census_design
from fetch_parent_tiles import tile_query


class StratumCensusTests(unittest.TestCase):
    def table(self, rows):
        return pd.DataFrame(rows,columns=['source_id','phot_g_mean_mag'])

    def test_census_replaces_sample_without_double_counting_and_preserves_boundary(self):
        faint=self.table([(1,12.),(2,18.99),(3,19.),(4,20.9)])
        bright=self.table([(5,11.99)])
        full=self.table([(1,12.),(2,18.99),(6,17.)])
        result=combine_parent_strata(faint,bright,full)
        self.assertEqual(set(result.source_id),{1,2,3,4,5,6})
        self.assertEqual(len(result),6)

    def test_invalid_range_and_duplicate_sources_fail(self):
        faint=self.table([(1,19.)]);bright=self.table([(2,11.)])
        with self.assertRaisesRegex(ValueError,'Supplement outside'):
            combine_parent_strata(faint,bright,self.table([(3,19.)]))
        with self.assertRaisesRegex(ValueError,'Repeated Gaia'):
            combine_parent_strata(faint,bright,self.table([(2,17.)]))

    def tiled(self):
        return dict(query_adql='indexed',g_min=12.,g_max=19.,random_sampling_modulus=1,
            radius_deg=1.,center=[10.,-20.],requested_tile_level=7,requested_tiles=[100],
            leaf_patches=[dict(pixels=[p],level=8,query_adql=tile_query([p],8,12.,19.,1))
                          for p in range(400,404)])

    def test_census_requires_exact_coverage_and_no_quality_cuts(self):
        model=SimpleNamespace(skymath=SimpleNamespace(cone_tiles=lambda *args:[100]))
        with patch('selection_bridge.load_study',return_value=model):
            verify_census_design(self.tiled(),10.,-20.,1.,12.,19.)
            bad=copy.deepcopy(self.tiled());bad['leaf_patches'].pop()
            with self.assertRaisesRegex(ValueError,'coverage'):
                verify_census_design(bad,10.,-20.,1.,12.,19.)
            bad=copy.deepcopy(self.tiled());bad['leaf_patches'][0]['query_adql']+=' AND ruwe<1.4'
            with self.assertRaisesRegex(ValueError,'selection design'):
                verify_census_design(bad,10.,-20.,1.,12.,19.)

    def test_original_two_stratum_design_is_unchanged(self):
        result=combine_parent_strata(self.table([(1,12.),(2,20.)]),self.table([(3,5.)]))
        self.assertEqual(result.source_id.tolist(),[1,2,3])


if __name__=='__main__':
    unittest.main()
