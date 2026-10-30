import json,tempfile,unittest
from pathlib import Path
import numpy as np
from population_diagnostics import population_summary,run,selected_arrays
from provenance import sha256_file
import test_joint_grid as fixtures
from joint_grid import analyse

class PopulationDiagnosticsTests(unittest.TestCase):
    def test_population_dependence_not_reported_as_confidence_or_point_measurement(self):
        rows=[dict(scenario=str(t),rate_msun_myr=r,duration_myr=1000.,negative_log_likelihood=(r-best)**2+offset,expected_signal_stars=100.)
              for t,best,offset in [(-1,31.,0.),(0,11.,.5),(1,5.,1.)] for r in [5.,11.,31.]]
        s=population_summary(rows)
        self.assertIsNone(s['primary_rate_msun_myr']);self.assertFalse(s['physical_rate_measurement_validated'])
        self.assertFalse(s['range_is_confidence_interval']);self.assertEqual(s['conditional_best_rate_range_msun_myr'],[5.,31.])
        self.assertAlmostEqual(s['conditional_best_rate_spread_factor'],6.2)
        by={r['scenario']:r for r in s['conditional_scenario_results']}
        self.assertEqual(by['1']['q_relative_to_global_min'],2.)
        self.assertAlmostEqual(by['1']['selected_stars_per_escaped_msun'],.02)
        self.assertTrue(by['-1']['rate_at_grid_boundary'])

    def test_native_grid_diagnostics_preserve_likelihood_and_detect_changed_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);paths=fixtures.JointGridTests().fixtures(root)
            analyse(paths,root/'scores',pixel_groups=np.array([0,0,1]))
            before=sha256_file(root/'scores/joint_grid.json')
            (root/'production_assessment.json').write_text(json.dumps(dict(best_sampled_candidate={'rate_msun_myr':2.},legacy_field='retained')))
            s=run(root,root/'diagnostics',compare_rates=(1.,2.))
            self.assertEqual(before,sha256_file(root/'scores/joint_grid.json'))
            assessment=json.loads((root/'production_assessment.json').read_text())
            self.assertIsNone(assessment['primary_rate_msun_myr']);self.assertEqual(assessment['legacy_field'],'retained')
            self.assertEqual(assessment['best_sampled_candidate']['rate_msun_myr'],2.)
            self.assertTrue((root/'diagnostics/previous_production_assessment.json').exists())
            self.assertFalse(s['changed_likelihood_weights']);self.assertFalse(s['external_population_constraint_applied'])
            with self.assertRaises(FileExistsError):run(root,root/'diagnostics')
            self.assertTrue((root/'diagnostics/population_visibility_grid.csv').exists())
            data=np.load(root/'diagnostics/candidate_0.npz')
            np.testing.assert_allclose(data['pearson_residual'],(data['observed']-data['expected'])/np.sqrt(data['expected']))
            self.assertEqual(len(data['observed']),2);data.close()
            score=json.loads((root/'scores/joint_grid.json').read_text());row=score['candidates'][0]
            Path(score['diagnostics'][0]['inputs'][0]['path']).write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError,'Scored template changed'):selected_arrays(score,row)

if __name__=='__main__':unittest.main()
