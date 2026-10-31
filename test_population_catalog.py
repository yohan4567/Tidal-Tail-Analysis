import unittest
from pathlib import Path
from population_catalog import build
from empirical_population import measured_constraint


class PopulationCatalogTests(unittest.TestCase):
    def test_real_source_registry_covers_sample_without_inventing_escape_constraints(self):
        root=Path(__file__).parent/'outputs/cluster_survey_29_v1'
        if not (root/'adopted_parameters.csv').exists():self.skipTest('Local survey snapshot not available')
        report=build(root);records=report['clusters']
        self.assertEqual(len(records),29)
        self.assertEqual(len({r['cluster'] for r in records}),29)
        byid={r['cluster']:r for r in records}
        self.assertEqual(byid['ESO452-SC11']['global_mf']['status'],'indirect_estimate_only')
        self.assertIsNone(byid['ESO452-SC11']['global_mf']['measured_mass_range_msun'])
        self.assertIsNone(byid['NGC5139']['age']['age_gyr'])
        self.assertTrue(byid['NGC6273']['age']['hold_reason'])
        self.assertTrue(all(not r['escaping_population']['total_absolute_rate_supported'] for r in records))
        self.assertTrue(all(r['age'].get('age_source') for r in records))

    def test_official_estimated_colour_is_rejected(self):
        html=Path(__file__).parent/'outputs/empirical_mf_source/parameter.html'
        if not html.exists():self.skipTest('Local official snapshot not available')
        with self.assertRaisesRegex(ValueError,'indirectly estimated'):
            measured_constraint(html,'ESO452-SC11')


if __name__=='__main__':unittest.main()
