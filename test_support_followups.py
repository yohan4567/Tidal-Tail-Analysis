import json,tempfile,unittest,zipfile
from pathlib import Path
from audit_grid_support import support_function
import portable_cluster
from collect_transfer_packages import collect
from provenance import sha256_file


def example_response(particles,draws):
    checked=[]
    for item in particles:
        checked.append(item)
        if item<0:raise ValueError('unsupported pair')
        for obs in range(draws):
            raise AssertionError('Observation loop must be omitted only in support audit')
    return checked


class SupportFollowupTests(unittest.TestCase):
    def test_support_audit_checks_every_batch_without_sampling_or_export(self):
        audit=support_function(example_response)
        self.assertEqual(audit([1,2,3],32),[1,2,3])
        with self.assertRaisesRegex(ValueError,'unsupported pair'):audit([1,-2,3],32)

    def test_portable_relocation_preserves_extended_error_library(self):
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp)/'bundle';base.mkdir();out=Path(temp)/'run'
            (base/'portable_manifest.json').write_text('{}')
            (base/'portable_cluster.py').write_text('fixture runner identity')
            original=base/'inputs/populations/offset_0';original.mkdir(parents=True)
            (original/'population_prior.json').write_text('{}')
            (original/'response_2.json').write_text(json.dumps(dict(errors={},completeness={})))
            info=dict(cluster='NGC6093',error_parent_file='inputs/calibration/parent_errors_extended.csv.gz',
                      error_parent_metadata='inputs/calibration/parent_errors_extended.csv.json')
            portable_cluster.prepare(base,out,info)
            rebound=json.loads((out/'offset_0/configs/response_2.json').read_text())
            self.assertEqual(Path(rebound['errors']['path']),base/'inputs/calibration/parent_errors_extended.csv.gz')
            self.assertEqual(Path(rebound['errors']['metadata']),base/'inputs/calibration/parent_errors_extended.csv.json')

    def test_transfer_copy_requires_validation_and_preserves_original(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);bundle=root/'NGC6304';bundle.mkdir();out=root/'transfer';out.mkdir()
            manifest=bundle/'portable_manifest.json';manifest.write_text('{}')
            archive=bundle.with_suffix('.zip')
            with zipfile.ZipFile(archive,'w') as z:z.writestr('NGC6304/portable_manifest.json','{}')
            validation=root/'validation.json';validation.write_text(json.dumps(dict(cluster='NGC6304',stage='relocated_small_response_passed',
                bundle_manifest_sha256=sha256_file(manifest))))
            entry=dict(archive=str(archive),validation=str(validation),analysis_kind='model_comparison',
                zip_crc_verified=True,archive_sha256=sha256_file(archive))
            result=collect(entry,'NGC6304',out)
            self.assertTrue(archive.is_file());self.assertEqual(sha256_file(result['archive']),sha256_file(archive))
            entry['zip_crc_verified']=False
            with self.assertRaisesRegex(ValueError,'validation'):collect(entry,'NGC6304',out)


if __name__=='__main__':unittest.main()
