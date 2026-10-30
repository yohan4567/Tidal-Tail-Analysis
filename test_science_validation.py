import json,tempfile,unittest
from pathlib import Path
from joint_grid import analyse
from science_validation import validate_grid,project_set

class ScienceValidationTests(unittest.TestCase):
    def test_production_binning_confidence_inversion(self):
        from test_joint_grid import JointGridTests
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);paths=JointGridTests().fixtures(root)
            analyse(paths,root/'scores',pixel_groups=[0,0,1])
            result=validate_grid(paths,root/'scores/joint_grid.json',root/'validation',20,20,123,[1.])
            self.assertEqual(result['confidence_level'],.95)
            self.assertEqual(len(result['coverage']),4)
            self.assertFalse(result['physical_rate_measurement_validated'])
            self.assertTrue((root/'validation/confidence_heatmap.png').exists())
            again=validate_grid(paths,root/'scores/joint_grid.json',root/'validation',20,20,123,[1.])
            self.assertEqual(result,again)
            with self.assertRaises(ValueError):validate_grid(paths,root/'scores/joint_grid.json',root/'validation',20,20,124,[1.])
    def test_disjoint_set_and_empty_set_not_fake_interval(self):
        result=project_set([1,2,3],[1,3])
        self.assertEqual(result['accepted_grid_values'],[1,3])
        self.assertEqual(result['envelope'],[1,3])
        self.assertTrue(result['touches_upper_grid_edge'])
        self.assertIsNone(project_set([1,2],[])['envelope'])

if __name__=='__main__':unittest.main()
