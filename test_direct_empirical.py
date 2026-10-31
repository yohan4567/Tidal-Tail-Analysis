import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_direct_empirical as module


class DirectQueueTests(unittest.TestCase):
    def test_preserves_grid_and_reuses_dynamics_without_legacy_grid(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);dest=root/'NGC2298';dest.mkdir()
            for name in ['calibration_inputs/isochrone.dat','science_calibration/bundle/tilt_plus0_zp_minus0.017.json',
                         'coupled_dynamics.json','selection/selection.json']:
                p=dest/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('{}')
            (root/'configs').mkdir();(root/'configs/survey_populations_v1.json').write_text(json.dumps(
                dict(clusters=dict(NGC2298=dict(age_gyr=12.67)))))
            def population(*args,**kw):
                out=args[4];out.mkdir(parents=True)
                return dict(response_configs=['fixture_response.json'])
            with patch.object(module,'HERE',root),patch.object(module,'measured_constraint'),\
                 patch.object(module,'prepare',side_effect=population),patch.object(module,'run_grid') as grid,\
                 patch.object(module,'infer') as infer:
                module.run(root,'NGC2298',12)
            self.assertEqual(grid.call_count,3)
            self.assertEqual(infer.call_count,3)
            calls=grid.call_args_list
            self.assertIsNone(calls[0].kwargs['reuse_grid'])
            self.assertEqual(calls[1].kwargs['reuse_grid'],calls[0].args[3])
            self.assertEqual(calls[2].kwargs['reuse_grid'],calls[0].args[3])
            for call in calls:
                self.assertEqual(call.kwargs['rates'],list(map(float,range(1,32,2))))
                self.assertEqual(call.kwargs['times'],list(map(float,range(100,2001,100))))
                self.assertEqual(len(call.kwargs['seeds']),3)
                self.assertEqual(call.kwargs['workers'],12)
                self.assertEqual(call.kwargs['response_workers'],12)
                self.assertEqual(call.kwargs['interval_myr'],.25)
            self.assertFalse((dest/'production_grid_shared_v1').exists())


if __name__=='__main__':unittest.main()
