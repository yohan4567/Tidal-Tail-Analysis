import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import run_deadline_queue as queue
import refresh_cluster_progress as progress


class Child:
    pid=321
    def __init__(self,code=0,sample=False):self.code=code;self.sample=sample
    def wait(self,timeout):
        if self.sample:
            self.sample=False
            raise subprocess.TimeoutExpired('fixture',timeout)
        return self.code


class ContinuousQueueTests(unittest.TestCase):
    def prepare(self,root,ids):
        (root/'configs').mkdir()
        (root/'configs/survey_populations_v1.json').write_text(json.dumps(
            {'clusters':{c:{'age_gyr':12.} for c in ids}}))
        for cid in ids:
            for name in queue.REQUIRED:
                p=root/cid/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('{}')

    def run_fixture(self,root,ids,children,**kw):
        with patch.object(queue,'HERE',root),patch.object(queue,'SAMPLE',[(c,c) for c in ids]),\
             patch.object(queue,'measured_constraint'),\
             patch.object(queue,'memory_snapshot',return_value=dict(total_gib=32.,available_gib=24.,used_gib=8.)),\
             patch.object(queue.subprocess,'Popen',side_effect=children) as popen:
            result=queue.run(root,12,**kw)
        return result,popen

    def test_child_error_continues_to_next_and_records_resources(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.prepare(root,['A','B'])
            state,popen=self.run_fixture(root,['A','B'],[Child(1),Child(sample=True)])
            self.assertEqual([e['stage'] for e in state['entries']],
                             ['blocked_or_failed','completed_conditional_posteriors_validation_pending'])
            self.assertEqual(state['stage'],'finished')
            self.assertEqual(popen.call_count,2)
            self.assertEqual(popen.call_args.args[0][-2:],['--workers','12'])
            events=[json.loads(line) for line in (root/'execution_events.jsonl').read_text().splitlines()]
            self.assertTrue(any(e['event']=='resource_sample' for e in events))
            self.assertTrue(any(e.get('exit_code')==1 for e in events))
            self.assertIn('blocked_or_failed',(root/'EXECUTION_LOG.md').read_text(encoding='utf-8'))

    def test_resume_prioritizes_interrupted_and_skips_failed_completed_and_pal5(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ids=['A','B','C','Pal5'];self.prepare(root,ids)
            (root/'deadline_queue_status.json').write_text(json.dumps({'entries':[
                {'cluster':'A','stage':'blocked_or_failed','reason':'known failure'},
                {'cluster':'B','stage':'running_direct_empirical'}]}))
            folder=root/'C/empirical_direct_v1';folder.mkdir()
            reports=[]
            for i in range(3):
                p=folder/f'report{i}.json';p.write_text('{}');reports.append(str(p))
            (folder/'transfer_sensitivity.json').write_text(json.dumps({'reports':reports}))
            state,popen=self.run_fixture(root,ids,[Child()])
            self.assertEqual(state['entries'][0]['cluster'],'B')
            stages={e['cluster']:e['stage'] for e in state['entries']}
            self.assertEqual(stages['A'],'blocked_or_failed')
            self.assertEqual(stages['C'],'skipped_completed')
            self.assertEqual(stages['Pal5'],'skipped_preserved')
            self.assertEqual(popen.call_count,1)

    def test_explicit_retry_attempts_failed_cluster_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.prepare(root,['A'])
            (root/'A/empirical_job.json').write_text('{"stage":"empirical_stage_failed"}')
            state,popen=self.run_fixture(root,['A'],[Child(1)],retry_failed=True)
            self.assertEqual(popen.call_count,1)
            self.assertEqual(state['entries'][0]['stage'],'blocked_or_failed')
            self.assertEqual(json.loads((root/'A/empirical_job.json').read_text())['stage'],'empirical_stage_failed')

    def test_ram_headroom_failure_does_not_start_child_and_does_not_stop_queue(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.prepare(root,['A','B'])
            with patch.object(queue,'HERE',root),patch.object(queue,'SAMPLE',[('A','A'),('B','B')]),\
                 patch.object(queue,'measured_constraint'),\
                 patch.object(queue,'memory_snapshot',return_value=dict(total_gib=32.,available_gib=1.,used_gib=31.)),\
                 patch.object(queue.subprocess,'Popen') as popen:
                state=queue.run(root,12,memory_budget_gib=24.)
            popen.assert_not_called()
            self.assertEqual(len(state['entries']),2)
            self.assertTrue(all(e['stage']=='blocked_or_failed' for e in state['entries']))

    def test_progress_prefers_queue_failure_over_stale_child_running_label(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'A').mkdir()
            (root/'A/empirical_job.json').write_text('{"stage":"running_empirical_response_grid"}')
            (root/'deadline_queue_status.json').write_text(json.dumps({'workers':12,'memory_budget_gib':24,
                'entries':[{'cluster':'A','stage':'blocked_or_failed','reason':'calibration coverage'}]}))
            with patch.object(progress,'SAMPLE',[('A','A')]):result=progress.refresh(root)
            self.assertEqual(result['entries'][0]['queue_job']['stage'],'blocked_or_failed')
            text=(root/'PROGRESS.md').read_text(encoding='utf-8')
            self.assertNotIn('running_empirical_response_grid',text)
            self.assertIn('calibration coverage',text)


if __name__=='__main__':unittest.main()
