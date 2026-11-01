import copy
import unittest
from types import SimpleNamespace
from unittest.mock import patch
from bright_calibration import verify_bright_design
from fetch_parent_tiles import tile_query
from selection_model import parent_query


class BrightTilesTests(unittest.TestCase):
    def test_original_cone_design_still_accepted(self):
        verify_bright_design({'query_adql':parent_query(10.,-20.,1.,0.,12.,1)},10.,-20.,1.)

    def tiled(self):
        leaves=[dict(pixels=[p],level=8,query_adql=tile_query([p],8,0.,12.,1)) for p in range(400,404)]
        return dict(query_adql='indexed',g_min=0.,g_max=12.,random_sampling_modulus=1,
                    radius_deg=1.,center=[10.,-20.],requested_tile_level=7,requested_tiles=[100],leaf_patches=leaves)

    def test_subdivided_disjoint_census_has_complete_coverage(self):
        model=SimpleNamespace(skymath=SimpleNamespace(cone_tiles=lambda *args:[100]))
        with patch('selection_bridge.load_study',return_value=model):
            verify_bright_design(self.tiled(),10.,-20.,1.)
            missing=self.tiled();missing['leaf_patches'].pop()
            with self.assertRaisesRegex(ValueError,'coverage'):verify_bright_design(missing,10.,-20.,1.)

    def test_no_quality_filtered_or_overlapping_census_accepted(self):
        model=SimpleNamespace(skymath=SimpleNamespace(cone_tiles=lambda *args:[100]))
        with patch('selection_bridge.load_study',return_value=model):
            bad=self.tiled();bad['leaf_patches'][0]['query_adql']+=' AND ruwe < 1.4'
            with self.assertRaisesRegex(ValueError,'selection design'):verify_bright_design(bad,10.,-20.,1.)
            bad=self.tiled();bad['leaf_patches'].append(copy.deepcopy(bad['leaf_patches'][0]))
            with self.assertRaisesRegex(ValueError,'Overlapping'):verify_bright_design(bad,10.,-20.,1.)


if __name__=='__main__':unittest.main()
