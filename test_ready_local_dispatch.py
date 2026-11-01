import json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from repair_seven_clusters import worker_lock
import dispatch_ready_local_repairs as dispatch

class ReadyDispatchTests(unittest.TestCase):
    def test_single_writer_lock(self):
        with tempfile.TemporaryDirectory() as d:
            with worker_lock(Path(d)):
                with self.assertRaises(OSError):
                    with worker_lock(Path(d)):self.fail('Second writer acquired lock')
            with worker_lock(Path(d)):pass
    def test_no_package_never_starts_numerical_work(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            with patch.object(dispatch,'preparing',return_value=False),patch.object(dispatch.subprocess,'run') as run:
                dispatch.run(root,12)
            run.assert_not_called()
            report=json.loads((root/'local_dispatch_status.json').read_text())
            self.assertTrue(all(r['stage']=='inputs_not_ready' for r in report['clusters'].values()))
    def test_zip_export_without_validation_is_not_dispatched(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            (root/'READY_PACKAGES.json').write_text(json.dumps({'ready_packages':{'NGC5904':{'zip_crc_verified':False}}}))
            with patch.object(dispatch,'preparing',return_value=False),patch.object(dispatch.subprocess,'run') as run:
                dispatch.run(root,12)
            run.assert_not_called()
            report=json.loads((root/'local_dispatch_status.json').read_text())
            self.assertEqual(report['clusters']['NGC5904']['stage'],'repaired_analysis_failed')

if __name__=='__main__':unittest.main()
