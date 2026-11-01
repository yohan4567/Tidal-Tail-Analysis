import json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from run_cluster_survey import SAMPLE
from run_deadline_queue import run,completed_reports

class DispatchTests(unittest.TestCase):
    def test_reserved_targets_never_start_a_child(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);ids=[c for _,c in SAMPLE]
            (root/'dispatch_overrides.json').write_text(json.dumps({'excluded_clusters':ids,'reason':'Other PC'}))
            with patch('run_deadline_queue.subprocess.Popen') as child:
                result=run(root,workers=1)
            child.assert_not_called()
            self.assertEqual({e['stage'] for e in result['entries']},{'reserved_for_other_pc'})
            self.assertEqual(len(result['entries']),29)

    def test_completion_requires_all_three_existing_reports(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);reports=[]
            for i in range(3):
                p=root/f'report_{i}.json';p.write_text('{}');reports.append(str(p))
            (root/'empirical_job.json').write_text(json.dumps({'stage':'completed_empirical_posteriors_validation_pending','reports':reports}))
            self.assertEqual(completed_reports(root),reports)
            Path(reports[-1]).unlink()
            self.assertEqual(completed_reports(root),[])

if __name__=='__main__':unittest.main()
