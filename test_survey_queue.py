import json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from provenance import sha256_file
from run_remaining_survey import checked_isochrone,run

class SurveyQueueTests(unittest.TestCase):
    def test_isochrone_reuse_requires_matching_age_abundance_and_checksum(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'iso.dat';self.assertFalse(checked_isochrone(p,1e10,-1.))
            p.write_text('official-input-test');m=dict(sha256=sha256_file(p),age_years=1e10,metallicity_mh=-1.)
            p.with_suffix('.json').write_text(json.dumps(m));self.assertTrue(checked_isochrone(p,1e10,-1.))
            with self.assertRaises(ValueError):checked_isochrone(p,2e10,-1.)
            p.write_text('changed')
            with self.assertRaises(ValueError):checked_isochrone(p,1e10,-1.)

    def test_cmd_clipped_model_is_rejected(self):
        import pandas as pd
        from fetch_isochrone import validate_returned_model
        table=pd.DataFrame(dict(MH=[-2.19174,-2.19174],logAge=[10.,10.]))
        with self.assertRaisesRegex(ValueError,'clipping/mismatch'):validate_returned_model(table,1e10,-2.31)
        validate_returned_model(table,1e10,-2.19174)
        with self.assertRaises(ValueError):validate_returned_model(table,1.1e10,-2.19174)

    def test_missing_population_is_reported_and_never_starts_a_grid(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'catalogs').mkdir();(root/'catalogs/master_clusters.csv').write_text('id\nNGC288\n')
            (root/'NGC288/selection').mkdir(parents=True);(root/'NGC288/selection/selection.json').write_text('{}')
            config=root/'ages.json';config.write_text('{"clusters":{}}')
            with patch('run_remaining_survey.DEFAULT_STUDY',root),patch('run_remaining_survey.execute') as execute:
                run(root,config,1,1,1,clusters=['NGC288']);execute.assert_not_called()
            status=json.loads((root/'NGC288/remaining_job.json').read_text())
            self.assertEqual(status['stage'],'blocked_population_review');self.assertFalse(status['physical_rate_measurement_validated'])
            self.assertFalse((root/'Pal5').exists())

if __name__=='__main__':unittest.main()
