import json,tempfile,unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import gammainc,gammaln
from likelihood_comparison import compare
from provenance import sha256_file


class LikelihoodComparisonTests(unittest.TestCase):
    def test_matches_poisson_integral_and_population_weighting_without_target_prior(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);grid=root/'grid';(grid/'scores').mkdir(parents=True);diagnostics=[]
            signals={('A',1.):.5,('B',1.):6.,('A',3.):3.,('B',3.):12.}
            for (scenario,rate),signal in signals.items():
                p=root/f'{scenario}_{rate}.npz';np.savez(p,counts=[5],background=[1.],signal_draws=np.array([[signal],[signal]]))
                diagnostics.append(dict(scenario=scenario,rate_msun_myr=rate,duration_myr=100.,inputs=[dict(path=str(p),sha256=sha256_file(p))]))
            (grid/'scores/joint_grid.json').write_text(json.dumps(dict(diagnostics=diagnostics,pixel_groups=None)))
            prior=root/'prior.json';prior.write_text(json.dumps(dict(scenario_weights={'A':.8,'B':.2})))
            report=compare(grid,root/'comparison',prior)
            f=pd.read_csv(root/'comparison/likelihood_grid.csv')
            for r in f.itertuples():
                exact=sum(weight*(gammainc(6,signals[(name,r.rate_msun_myr)]+5)-gammainc(6,signals[(name,r.rate_msun_myr)]))/5 for name,weight in [('A',.8),('B',.2)])
                self.assertAlmostEqual(r.log_likelihood,np.log(exact)+gammaln(6),places=7)
            self.assertFalse(report['target_parameter_posterior_computed'])
            self.assertFalse(report['confidence_intervals_calibrated'])
            self.assertEqual(report['grid_points'],2)


if __name__=='__main__':unittest.main()
