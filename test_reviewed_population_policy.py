import unittest
from pathlib import Path
from repair_four_populations import global_metallicity,builder,OLD
from reviewed_population_policy import prepare,indirect_constraint

class ReviewedPopulationTests(unittest.TestCase):
    def test_no_alpha_enhancement_keeps_feh(self):
        self.assertAlmostEqual(global_metallicity(-2.35,0),-2.35)
        self.assertGreater(global_metallicity(-2.35,.46),-2.2)
    def test_indirect_mf_is_not_relabelled_as_measured(self):
        fixed=builder()
        with self.assertRaises(ValueError):fixed.measured_constraint(OLD/'mf_source.html','ESO452-SC11')
        c=indirect_constraint(OLD/'mf_source.html','ESO452-SC11')
        self.assertTrue(c['indirect_estimate']);self.assertIsNone(c['measured_mass_range_msun'])
    def test_eso_still_requires_explicit_opt_in(self):
        with self.assertRaisesRegex(ValueError,'opt-in'):prepare('ESO452-SC11',None,None,None,None)
        with self.assertRaisesRegex(ValueError,'scoped'):prepare('NGC6341',None,None,None,None,allow_indirect_model=True)

if __name__=='__main__':unittest.main()
