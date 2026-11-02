import unittest
from repair_sparse_strata import progress_record


class PreparationProgressTests(unittest.TestCase):
    def test_verified_package_entry_can_be_reported_ready(self):
        entry=dict(archive='NGC6355_v3.zip',zip_crc_verified=True,
                   physical_rate_measurement_validated=False)
        result=progress_record({},'conditional_model_package_ready',**entry)
        self.assertEqual(result['stage'],'conditional_model_package_ready')
        self.assertEqual(result['archive'],entry['archive'])
        self.assertTrue(result['zip_crc_verified'])
        self.assertFalse(result['physical_rate_measurement_validated'])

    def test_recovered_status_clears_old_failure(self):
        result=progress_record(dict(error='old status-write failure'),'conditional_model_package_ready')
        self.assertNotIn('error',result)

    def test_ready_package_does_not_imply_physical_rate_validation(self):
        result=progress_record({},'conditional_model_package_ready',physical_rate_measurement_validated=True)
        self.assertFalse(result['physical_rate_measurement_validated'])


if __name__=='__main__':unittest.main()
