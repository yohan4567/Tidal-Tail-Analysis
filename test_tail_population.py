import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from population_model import Population
from provenance import sha256_file
from tail_population import prepare, read_constraint


class TailPopulationTests(unittest.TestCase):
    def fixture(self, root):
        constraint = root/'constraint.json'
        constraint.write_text(json.dumps(dict(cluster='fixture', source='test fixture', citation='fixture',
            convention='dN/dm proportional to m**(-(1+x))', mass_coordinate='main_sequence_current_mass',
            measurement_scope='present_day_observed_tail_segment', mass_range_msun=[.58,.82],
            regions=dict(south=dict(x=.97,sigma_x=.36)))))
        iso = root/'iso.csv'
        pd.DataFrame(dict(Mini=[.1,.2,.4,.6,.8,.81,.82], Mass=[.1,.2,.4,.6,.8,.7,.6],
            int_IMF=np.arange(7.), Gmag=np.arange(7.), G_BPmag=np.arange(7.)+.3,
            G_RPmag=np.arange(7.)-.5, label=[1,1,1,1,1,3,4])).to_csv(iso,index=False)
        iso.with_suffix('.json').write_text(json.dumps(dict(sha256=sha256_file(iso),age_years=1e10,metallicity_mh=-1.)))
        ref=root/'response.json'
        ref.write_text(json.dumps(dict(population='old', quality='quality.json', cache_dir='cache',
            errors=dict(path='errors.csv',metadata='errors.json'), completeness=dict(path='complete.csv',metadata='complete.json'),
            parallax_zero_point_mas=-.017)))
        return constraint,iso,ref

    def test_literature_sign_conversion(self):
        path=Path(__file__).parent/'configs/Pal5_tail_mf_Koch2004.json'
        self.assertAlmostEqual(read_constraint(path,'south')['alpha_mean'],-1.97)
        self.assertAlmostEqual(read_constraint(path,'north')['alpha_mean'],-2.31)

    def test_no_silent_extrapolation_or_partial_output(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);c,i,r=self.fixture(root)
            with self.assertRaisesRegex(ValueError,'does not cover'):
                prepare(c,'south',i,r,root/'out')
            self.assertFalse((root/'out').exists())

    def test_conditional_population_moments_mass_and_provenance(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);c,i,r=self.fixture(root)
            report=prepare(c,'south',i,r,root/'out',allow_extrapolation=True)
            self.assertFalse(report['total_absolute_rate_supported'])
            scenarios=[json.loads(s) for s in report['scenario_weights']]
            w=np.array(list(report['scenario_weights'].values()))
            a=np.array([s['measured_tail_mf_slope'] for s in scenarios])
            self.assertAlmostEqual(w.sum(),1.)
            self.assertAlmostEqual(w@a,-1.97)
            self.assertAlmostEqual(w@((a+1.97)**2),.36**2)
            self.assertTrue(all('escape_slope_offset' not in s for s in scenarios))
            pop=Population.load(root/'out/population_0')
            self.assertAlmostEqual(float((pop.table.current_mass_msun*pop.table.number_per_msun).sum()),1.)
            self.assertTrue(pop.metadata['support_audit']['mf_extrapolated'])
            self.assertEqual(report['constraint']['constraint_sha256'],sha256_file(c))
            with self.assertRaisesRegex(ValueError,'fresh'):
                prepare(c,'south',i,r,root/'out',allow_extrapolation=True)

    def test_refuse_unknown_sign_convention_and_changed_isochrone(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);c,i,r=self.fixture(root)
            i.write_text(i.read_text()+'\n')
            with self.assertRaisesRegex(ValueError,'checksum'):
                prepare(c,'south',i,r,root/'out',allow_extrapolation=True)
            x=json.loads(c.read_text());x['convention']='Salpeter';c.write_text(json.dumps(x))
            with self.assertRaisesRegex(ValueError,'convention'):
                read_constraint(c,'south')


if __name__=='__main__':
    unittest.main()
