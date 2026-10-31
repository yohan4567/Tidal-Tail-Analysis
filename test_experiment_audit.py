import json
from pathlib import Path
import tempfile
import unittest

from audit_experiment import audit_posterior


class ExperimentAuditTests(unittest.TestCase):
    def test_measured_tail_and_narrow_intervals_do_not_imply_total_mass_closure(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'posterior.json'
            axis=dict(median=5.,credible68=[4.9,5.1],credible95=[4.8,5.2],
                first_node_quadrature_mass=0.,last_node_quadrature_mass=0.)
            p.write_text(json.dumps(dict(population_prior=dict(
                constraint=dict(measurement_scope='present_day_observed_tail_segment'),
                support_audits=[dict(mf_extrapolated=True)]),
                results=dict(a=dict(rate=axis,lookback_time=axis,support=dict(rate_msun_myr=[1,31],time_myr=[100,2000]))))))
            result=audit_posterior(p)
            self.assertTrue(result['measured_tail_constraint'])
            self.assertTrue(result['mf_extrapolation_detected'])
            self.assertFalse(result['total_absolute_rate_supported'])
            self.assertIsNone(result['primary_total_absolute_rate_msun_myr'])
            self.assertEqual(next(c for c in result['checks'] if c['stage']=='zero_signal_model')['status'],'not_computed')


if __name__=='__main__':unittest.main()
