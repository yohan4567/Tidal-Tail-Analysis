import tempfile
import unittest
from pathlib import Path
from plot_mlr_parameters import comparison_directories


class ResultDiscoveryTests(unittest.TestCase):
    def test_new_returned_clusters_are_found_at_each_offset(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for name in ['12345','NGC 12346']:
                for offset in ['0','-0.5','-1']:
                    p=root/name/offset/'comparison';p.mkdir(parents=True)
                    (p/'comparison.json').write_text('{}')
            for d in [0.,-.5,-1.]:
                found=comparison_directories(d,root)
                self.assertEqual(found['NGC12345'],root/'12345'/f'{d:g}'/'comparison')
                self.assertIn('NGC12346',found)

    def test_duplicate_versions_are_not_silently_selected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for sub in ['0','offset_0']:
                p=root/'12345'/sub/'comparison';p.mkdir(parents=True)
                (p/'comparison.json').write_text('{}')
            with self.assertRaisesRegex(ValueError,'Multiple returned'):
                comparison_directories(0.,root)


if __name__=='__main__':unittest.main()
