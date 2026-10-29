import unittest
import json
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from dense_convergence import thinning_multipliers,tv


class DenseTests(unittest.TestCase):
    def test_full_analysis_builds_all_density_levels_before_comparison(self):
        from dense_convergence import analyse
        from response_windows import rebin
        from provenance import sha256_file
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);study=root/'study';selection=root/'selection'
            study.mkdir();selection.mkdir()
            (study/'plan.json').write_text(json.dumps(dict(rates=[9.],seeds=[1,2,3],release_interval_myr=.25,times_myr=[250,500,750,1000])))
            np.savez(selection/'maps.npz',counts=np.ones((4,4)),background=np.ones((4,4)),fitmask=np.ones((4,4),bool),edges=np.linspace(-4,4,5))
            ages=np.tile((np.arange(4000)+.5)*.25,2);pixels=np.arange(8000)%16
            contributions=np.full((2,8000),.001)
            for seed in (1,2,3):
                parent=study/f'rate_9_seed_{seed}';(parent/'simulation').mkdir(parents=True)
                pd.DataFrame(dict(lookback_release_myr=ages)).to_csv(parent/'simulation/particles.csv',index=False)
                full=parent/'population_0_T_1000';full.mkdir()
                np.savez(full/'calibrated_counts.npz',signal_draws=rebin(contributions,pixels,ages,1000,16))
                (full/'response.json').write_text(json.dumps(dict(template_sha256=sha256_file(full/'calibrated_counts.npz'))))
                np.savez(full/'tracer_response.npz',lookback_myr=ages,calibration_counts=contributions,fit_pixel_index=pixels,mass_weight_msun=np.full(8000,1.125))
            result=analyse(study,selection,root/'analysis')
            self.assertEqual(len(result['noise']),20)
            self.assertTrue(all(r['passed'] for r in result['normalization']))
            self.assertEqual(len(pd.read_csv(root/'analysis/nested_thinning.csv')),40)

    def test_nested_stratified_weights_conserve_each_arm_block(self):
        levels=thinning_multipliers(8000,123)
        for density,m in levels.items():
            self.assertEqual(np.count_nonzero(m),density)
            np.testing.assert_array_equal(m.reshape(-1,4).sum(axis=1),4)
        self.assertTrue(np.all((levels[2000]>0)<=(levels[4000]>0)))
        with self.assertRaises(ValueError):thinning_multipliers(6,1)

    def test_tv_is_shape_only(self):
        self.assertEqual(tv([1,2],[2,4]),0)
        self.assertEqual(tv([1,0],[0,1]),1)


if __name__=='__main__':unittest.main()
