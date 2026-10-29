import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from joint_grid import window_weights, fixed_rate_fit, noise_report, analyse, read_templates
from joint_recovery import recover


class JointGridTests(unittest.TestCase):
    def test_multimodal_background_profile_finds_global_peak(self):
        from joint_fit import fixed_rate_fit as global_fit
        n=[14,9,10]; b=[.20106681,.4614672,2.09966255]
        s=[[14.09604319,.05953345,.1734261],[24.72208606,12.40441953,5.78438578],[3.82114504,2.01926236,5.94827532]]
        beta,value=global_fit(n,s,b)
        self.assertLess(value,-41.97)
        self.assertAlmostEqual(beta,1.5324,places=3)

    def test_rebin_windows_match_explicit_particle_sums(self):
        from response_windows import rebin
        states=np.array([[1.,2.,3.,4.],[2.,3.,4.,5.]])
        pixels=np.array([0,1,0,-1]); ages=np.array([.5,1.5,2.5,.5])
        np.testing.assert_array_equal(rebin(states,pixels,ages,2,2),[[1,2],[2,3]])
        np.testing.assert_array_equal(rebin(states,pixels,ages,3,2),[[4,2],[6,3]])

    def test_binomial_limits_remain_nonzero_at_boundary(self):
        from joint_recovery import wilson
        self.assertGreater(wilson(0,100)[1],0)
        self.assertLess(wilson(100,100)[0],1)

    def test_subset_keeps_mass_and_rejects_partial_stratum(self):
        particles=pd.DataFrame(dict(lookback_release_myr=[.2,1.2,2.2,3.2,.7,1.7,2.7,3.7]))
        run=dict(config=dict(duration_myr=4,n_particles=8,reference_mdot_msun_myr=3),
            release_time_sampling='stratified_lookback_both_arms',dynamics_mass_history=dict(mode='prescribed_linear_escape'))
        w,meta=window_weights(particles,run,2,np.full(8,1.5))
        self.assertEqual(w.sum(),6); self.assertEqual(meta['expected_mass_msun'],6)
        np.testing.assert_array_equal(w,[1.5,1.5,0,0,1.5,1.5,0,0])
        with self.assertRaises(ValueError): window_weights(particles,run,2.5,np.ones(8))
        run['dynamics_mass_history']=None
        with self.assertRaises(ValueError): window_weights(particles,run,2,np.ones(8))

    def test_background_only_and_zero_data_supported(self):
        beta,_=fixed_rate_fit([2,4],[[0,0]],[1,2]); self.assertAlmostEqual(beta,2,places=6)
        beta,_=fixed_rate_fit([0,0],[[0,0]],[1,2]); self.assertEqual(beta,0)
        with self.assertRaises(ValueError): fixed_rate_fit([1.5],[[1]],[1])

    def test_remote_region_preserves_covariance(self):
        reps=np.array([[1.,1.],[3.,3.],[5.,5.]])
        r=noise_report(reps,np.array([100.,100.]),{'remote':np.array([True,True])})
        expected=reps.sum(axis=1).std(ddof=1)/np.sqrt(3)
        self.assertAlmostEqual(r['regions']['remote']['standard_error'],expected)
        self.assertEqual(noise_report(reps[:1],np.array([100.,100.]))['status'],'unverified')

    def fixtures(self,root):
        paths=[]
        for duration in (1.,2.):
            for rate in (1.,2.):
                for seed in (42,84):
                    path=root/f'T{duration}_R{rate}_S{seed}.npz'
                    desc=dict(observation_identity='obs',selection_sha256='selection',mass_basis='living',
                        model_version='prescribed_mass_history_v1',dynamics_family_identity='family',
                        amplitude_policy='fixed_reference_rate',dynamics_mass_history=dict(stellar_rate_msun_myr=rate),
                        lookback_window=dict(duration_myr=duration),scenario=dict(escape_mass_tilt=0),
                        response_sampling_identity='aligned',dynamics_sampling_identity=str(seed))
                    # Saturation: the older bin contributes no stars to the selected footprint.
                    np.savez(path,counts=[12,4,2],background=[1.,1.,1.],signal_draws=np.array([[5.,1.,0.]])*rate,
                        reference_mdot=rate,selection_description=json.dumps(desc))
                    paths.append(path)
        return paths

    def test_joint_profile_plateau_does_not_claim_lower_limit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); paths=self.fixtures(root)
            report=analyse(paths,root/'scores')
            profile=pd.read_csv(root/'scores/time_profile.csv')
            np.testing.assert_allclose(profile.q,0,atol=1e-9)
            self.assertFalse(report['physical_rate_measurement_validated'])
            self.assertEqual(len(report['candidates']),4)
            with self.assertRaises(ValueError): analyse([paths[0],paths[0]],root/'bad')

    def test_recovery_held_out_and_null(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); paths=self.fixtures(root)
            manifest=root/'injections.json'
            manifest.write_text(json.dumps(dict(source_description='Software regression only',scenarios=[
                dict(name='interior',template=paths[0].name,tags=['interior']),
                dict(name='background',null=True,tags=['null'])])))
            report=recover(paths,manifest,root/'recovery',trials=40)
            self.assertEqual(report['scenarios'][0]['held_out_trials'],20)
            self.assertIsNotNone(report['scenarios'][0]['coverage'])
            self.assertIsNone(report['scenarios'][1]['truth_key'])
            self.assertIn('off_grid',report['missing_required_scenario_tags'])
            self.assertFalse(report['physical_measurement_validated'])


if __name__=='__main__': unittest.main()
