import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.special import xlogy

from fit_diagnostics import diagnose,mixture_deviance
from validate_science import absolute_spec,numerical_comparison,clear_empty_interrupted_simulation,format_change,check_response
from provenance import sha256_file


class ValidationChecks(unittest.TestCase):
    def test_saturated_deviance_includes_zero_counts(self):
        n=np.array([0.,3.,5.]); nll=float(np.sum(n-xlogy(n,n)))
        self.assertAlmostEqual(mixture_deviance(n,nll),0.)
        self.assertAlmostEqual(mixture_deviance(n,nll+2),4.)

    def test_refitted_bootstrap_reproducible_and_bounded(self):
        args=([20,8,6],[[3,1,0],[2.8,1.2,0]],[1,1,1])
        result,rows,pred=diagnose(*args,trials=20,seed=7)
        again,rows2,_=diagnose(*args,trials=20,seed=7)
        np.testing.assert_array_equal(rows.deviance,rows2.deviance)
        self.assertEqual(result['p_conditional_plus_one'],again['p_conditional_plus_one'])
        self.assertGreaterEqual(result['p_conditional_plus_one'],1/21)
        self.assertLessEqual(result['p_conditional_plus_one'],1)
        self.assertEqual(len(rows),20); self.assertTrue((pred>0).all())
        self.assertTrue(set(rows.calibration_state).issubset({0,1}))

    def test_cloned_response_config_keeps_absolute_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); p=root/'config.json'
            p.write_text(json.dumps(dict(population='population',quality='quality',cache_dir='../cache',
                completeness=dict(path='map.fits'),errors=dict(path='parent.csv.gz',metadata='parent.json'))))
            spec=absolute_spec(p)
            self.assertEqual(spec['population'],str((root/'population').resolve()))
            self.assertEqual(spec['errors']['path'],str((root/'parent.csv.gz').resolve()))

    def test_numerical_check_tracks_rate_and_spatial_shape(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for name,s in [('a',[1,1,0]),('b',[1,0,1])]:
                (root/name).mkdir(); np.savez(root/name/'calibrated_counts.npz',counts=[3,2,1],signal=s)
            a=dict(name='a',response=str(root/'a'),rate=dict(mdot_msun_myr=10.,reference_selected_stars=2.))
            b=dict(name='b',response=str(root/'b'),rate=dict(mdot_msun_myr=11.,reference_selected_stars=2.))
            check=numerical_comparison(a,b,.05)
            self.assertFalse(check['within_rate_and_count_tolerance'])
            self.assertAlmostEqual(check['normalized_signal_total_variation'],.5)

    def test_empty_interrupted_output_recovery_preserves_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); dest=root/'simulation'; dest.mkdir()
            clear_empty_interrupted_simulation(dest,root)
            self.assertFalse(dest.exists())
            dest.mkdir(); (dest/'particles.csv').write_text('partial data')
            with self.assertRaises(ValueError): clear_empty_interrupted_simulation(dest,root)
            self.assertEqual((dest/'particles.csv').read_text(),'partial data')
            self.assertIn('undefined',format_change(None))

    def test_cached_completeness_map_must_still_match(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); survey=root/'map.fits'; survey.write_bytes(b'map-version-1')
            template=root/'calibrated_counts.npz'; template.write_bytes(b'template')
            meta=dict(input_sha256={},source_sha256={},versions={},template_sha256=sha256_file(template),
                completeness=dict(path=str(survey),input_sha256=sha256_file(survey)))
            (root/'response.json').write_text(json.dumps(meta))
            check_response(root)
            survey.write_bytes(b'map-version-2')
            with self.assertRaisesRegex(ValueError,'Completeness map changed'): check_response(root)


if __name__=='__main__': unittest.main()
