import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import queue_repaired_seven as queue
from prepare_exact_local_retries import physical_config
from provenance import sha256_file


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding='utf-8')


class SevenQueueTests(unittest.TestCase):
    def test_waits_for_both_local_dispatch_and_ngc3201(self):
        with tempfile.TemporaryDirectory() as temp:
            compact, bright = Path(temp) / 'compact', Path(temp) / 'bright'
            self.assertFalse(queue.predecessors_finished(compact, bright))
            save(compact / 'local_dispatch_status.json', {'clusters': {
                cid: {'stage': queue.COMPLETE} for cid in queue.TARGETS[:2]}})
            save(bright / 'status.json', {'stage': 'running_repaired_model_comparisons'})
            self.assertFalse(queue.predecessors_finished(compact, bright))
            save(bright / 'status.json', {'stage': queue.COMPLETE})
            self.assertTrue(queue.predecessors_finished(compact, bright))

    def package(self, root, cid):
        archive = root / (cid + '.zip')
        archive.parent.mkdir(parents=True, exist_ok=True)
        archive.write_bytes(b'already CRC checked archive fixture')
        manifest = archive.with_suffix('') / 'portable_manifest.json'
        save(manifest, {'cluster': cid, 'analysis_kind': 'model_comparison'})
        validation = root / (cid + '_validation.json')
        save(validation, {'stage': 'relocated_small_response_passed',
                          'bundle_manifest_sha256': sha256_file(manifest)})
        return dict(archive=str(archive), archive_sha256=sha256_file(archive),
                    validation=str(validation), zip_crc_verified=True,
                    analysis_kind='model_comparison', job=str(root / (cid + '_job.json')))

    def test_rejects_changed_archive_or_relocation_identity(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            entry = self.package(root, 'NGC6341')
            queue.checked_entry(entry, 'NGC6341')
            Path(entry['archive']).write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'changed'):
                queue.checked_entry(entry, 'NGC6341')
            entry['archive_sha256'] = sha256_file(entry['archive'])
            save(Path(entry['archive']).with_suffix('') / 'portable_manifest.json', {'changed': True})
            with self.assertRaisesRegex(ValueError, 'manifest'):
                queue.checked_entry(entry, 'NGC6341')

    def test_failure_continues_and_completed_predecessors_are_not_rerun(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            root, compact, network, bright = (base / name for name in ('queue', 'compact', 'network', 'bright'))
            root.mkdir()
            previous = {}
            for cid in queue.TARGETS[:2]:
                final = compact / cid / 'run_status.json'
                save(final, {'stage': queue.COMPLETE})
                previous[cid] = {'stage': queue.COMPLETE, 'result': str(final)}
            save(compact / 'local_dispatch_status.json', {'clusters': previous})
            save(bright / 'status.json', {'stage': queue.COMPLETE})
            entries = {cid: self.package(compact, cid) for cid in ('NGC6341', 'NGC7099')}
            save(compact / 'READY_PACKAGES.json', {'ready_packages': entries})
            save(compact / 'finished.json', {})
            save(network / 'finished.json', {})
            def execute(args, **kwargs):
                job = args[args.index('--job') + 1]
                if 'NGC6341' in job:
                    return SimpleNamespace(returncode=1)
                out = Path(args[args.index('--out') + 1])
                save(out / 'run_status.json', {'stage': queue.COMPLETE})
                return SimpleNamespace(returncode=0)
            with patch.object(queue.subprocess, 'run', side_effect=execute) as process, \
                    patch('run_deadline_queue.memory_snapshot', return_value={'available_gib': 16}):
                queue.run(root, compact, network, bright)
            self.assertEqual(process.call_count, 2)
            states = queue.read(root / 'status.json')['clusters']
            self.assertEqual(states['NGC6341']['stage'], 'analysis_failed')
            self.assertEqual(states['NGC7099']['stage'], queue.COMPLETE)
            self.assertTrue(states['NGC5904']['reused_completed_run'])
            self.assertEqual(states['NGC6266']['stage'], 'inputs_failed')
            with patch.object(queue.subprocess, 'run') as repeated:
                queue.run(root, compact, network, bright)
            repeated.assert_not_called()

    def test_exact_reuse_ignores_only_annotation_not_physics(self):
        a = {'mass_evolution': {'source_description': 'original', 'rate_msun_myr': 5}, 'n_particles': 16000}
        b = {'mass_evolution': {'source_description': 'new description', 'rate_msun_myr': 5}, 'n_particles': 16000}
        self.assertEqual(physical_config(a), physical_config(b))
        b['mass_evolution']['rate_msun_myr'] = 7
        self.assertNotEqual(physical_config(a), physical_config(b))
        self.assertIn('source_description', a['mass_evolution'])


if __name__ == '__main__':
    unittest.main()
