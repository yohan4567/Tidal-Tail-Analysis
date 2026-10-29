import unittest
import json
import tempfile
from pathlib import Path
import numpy as np
from profile_recovery import rank_threshold,statistics


class ProfileRecoveryTests(unittest.TestCase):
    def test_small_frozen_template_recovery_keeps_targets_separate(self):
        from test_joint_grid import JointGridTests
        from profile_recovery import recover
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);paths=JointGridTests().fixtures(root)
            manifest=root/'injections.json'
            manifest.write_text(json.dumps(dict(source_description='Software test only',scenarios=[
                dict(name='signal',duration_myr=1,rate_msun_myr=1,scenario=dict(escape_mass_tilt=0),inject_exact_grid_mean=True),
                dict(name='null',null=True)])))
            result=recover(paths,manifest,root/'result',calibration=20,validation=20)
            targets=result['scenarios'][0]['statistics']
            self.assertEqual(set(targets),{'rate_profile','time_profile','time_rate_joint','candidate_point'})
            self.assertEqual(targets['rate_profile']['calibration_rank'],20)
            self.assertEqual(set(result['scenarios'][1]['statistics']),{'null_search'})
            self.assertFalse(result['physical_rate_measurement_validated'])

    def test_rank_correction_and_small_sample(self):
        threshold,k=rank_threshold(np.arange(50))
        self.assertEqual(k,49);self.assertEqual(threshold,48)
        self.assertEqual(np.quantile(np.arange(50),.95,method='higher'),47)
        self.assertTrue(np.isinf(rank_threshold(np.arange(10))[0]))
        with self.assertRaises(ValueError):rank_threshold([np.nan])

    def test_nuisance_parameters_profiled_not_fixed(self):
        keys=[(1,4,'a'),(1,4,'b'),(2,4,'a'),(1,9,'a')]
        q=statistics([6.,3.,2.,1.],keys,(1,4,'a'))
        self.assertEqual(q['candidate_point'],10.)
        self.assertEqual(q['time_rate_joint'],4.)
        self.assertEqual(q['time_profile'],0.)
        self.assertEqual(q['rate_profile'],2.)


if __name__=='__main__':unittest.main()
