import json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from transfer_status_io import write_json

class StatusTests(unittest.TestCase):
    def test_transient_replace_lock_retries(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'status.json';original=Path.replace;calls=[]
            def replace(source,target):
                calls.append(1)
                if len(calls)<3:raise PermissionError('Sharing violation')
                return original(source,target)
            with patch.object(Path,'replace',replace):write_json(p,{'stage':'ready'},delay=0)
            self.assertEqual(json.loads(p.read_text()),{'stage':'ready'})
            self.assertEqual(len(calls),3)
    def test_progress_fallback_does_not_apply_to_scientific_manifest(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'status.json'
            with patch.object(Path,'replace',side_effect=PermissionError('Locked')):
                write_json(p,{'stage':'ready'},status_file=True,attempts=1)
                with self.assertRaises(PermissionError):write_json(Path(d)/'manifest.json',{'value':1},attempts=1)
            self.assertEqual(json.loads(p.read_text())['stage'],'ready')

if __name__=='__main__':unittest.main()
