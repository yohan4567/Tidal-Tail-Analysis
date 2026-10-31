import json,tempfile,unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import gammaln,gammainc
from bayesian_grid import marginal_background,posterior_surface,density_quantile,integration_weights
from empirical_population import measured_constraint,normal_nodes,measured_population

class BayesianInferenceTests(unittest.TestCase):
    def test_background_only_integral_far_from_prior_support_does_not_underflow(self):
        from scipy.integrate import quad
        upper=.01;n=100.
        result=marginal_background([n],[[0.]],[1.],upper)
        reference=n*np.log(upper)+np.log(quad(lambda u:u**n*np.exp(-upper*u),0,1)[0])
        self.assertAlmostEqual(result['log_marginal_likelihood'],reference,places=8)

    def test_background_integral_matches_truncated_gamma_analytic_case(self):
        n=np.array([2.,3.]);b=np.array([1.,2.]);upper=5.;r=marginal_background(n,[[0.,0.]],b,upper)
        shape=n.sum()+1;B=b.sum()
        expected=(n*np.log(b)).sum()+gammaln(shape)-shape*np.log(B)+np.log(gammainc(shape,B*upper))-np.log(upper)
        self.assertAlmostEqual(r['log_marginal_likelihood'],expected,places=8)
        mean=shape/B*gammainc(shape+1,B*upper)/gammainc(shape,B*upper)
        self.assertAlmostEqual(r['background_mean'],mean,places=7)
    def test_continuous_quantiles_and_nonuniform_grid_mass(self):
        x=np.array([1.,2.,5.]);self.assertAlmostEqual(integration_weights(x).sum(),4.)
        self.assertAlmostEqual(density_quantile(x,np.ones(3),.16),1.64)
        nodes=[dict(scenario='a',rate_msun_myr=r,duration_myr=t,log_marginal_likelihood=0.) for r in x for t in [100.,200.,500.]]
        summary,a=posterior_surface(nodes,{'a':1.})
        self.assertAlmostEqual(summary['rate']['median'],3.)
        self.assertAlmostEqual(summary['log_evidence_up_to_common_count_factorial'],0.)
        for rp,tp in [('uniform','log_uniform'),('log_uniform','uniform'),('log_uniform','log_uniform')]:
            s,_=posterior_surface(nodes,{'a':1.},rp,tp)
            self.assertAlmostEqual(s['log_evidence_up_to_common_count_factorial'],0.)
        self.assertAlmostEqual(np.sum(a['joint_density']*integration_weights(x)[:,None]*integration_weights(a['times'])[None,:]),1.)
        nodes.pop()
        with self.assertRaises(ValueError):posterior_surface(nodes,{'a':1.})
    def test_population_prior_is_applied_before_marginalization(self):
        nodes=[dict(scenario=k,rate_msun_myr=r,duration_myr=t,log_marginal_likelihood=(-r if k=='low' else r)) for k in ['low','high'] for r in [1.,2.,3.] for t in [100.,200.]]
        low,_=posterior_surface(nodes,{'low':.99,'high':.01});high,_=posterior_surface(nodes,{'low':.01,'high':.99})
        self.assertLess(low['rate']['median'],high['rate']['median'])
        self.assertIn('credible95',low['rate'])
    def test_end_to_end_native_templates_and_completed_resume(self):
        import test_joint_grid as fixtures
        from joint_grid import analyse
        from bayesian_grid import run
        from provenance import sha256_file
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);paths=fixtures.JointGridTests().fixtures(root)
            analyse(paths,root/'scores',pixel_groups=np.array([0,0,1]))
            before=sha256_file(root/'scores/joint_grid.json')
            report=run(root,root/'posterior')
            self.assertEqual(sha256_file(root/'scores/joint_grid.json'),before)
            self.assertFalse(report['population_combined']);self.assertEqual(len(report['results']),4)
            self.assertEqual(report['background_only_model']['status'],'computed')
            self.assertTrue(all(np.isfinite(v['log_bayes_factor_vs_background_only']) for v in report['results'].values()))
            self.assertFalse(report['credible_intervals_are_confidence_intervals'])
            self.assertTrue((root/'posterior/joint_posterior.png').exists())
            again=run(root,root/'posterior');self.assertEqual(again['results'],report['results'])
            with self.assertRaises(ValueError):run(root,root/'posterior',beta_upper=2.)

    def test_mf_quadrature_and_sign(self):
        x,w=normal_nodes(-.84,.19,5);self.assertAlmostEqual(w.sum(),1.)
        self.assertAlmostEqual(np.sum(x*w),-.84);self.assertAlmostEqual(np.sum(w*(x+.84)**2),.19**2)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'table.html';cells=['Pal 5']+['x']*18+['0.29 - 0.81','-0.84 +/- 0.19']
            p.write_text('<table><tr>'+''.join('<td>'+v+'</td>' for v in cells)+'</tr></table>')
            m=measured_constraint(p,'Pal5');self.assertEqual(m['slope_mean'],-.84)
            p.write_text(p.read_text().replace('<td>-0.84','<td style="color:green">-0.84'))
            with self.assertRaises(ValueError):measured_constraint(p,'Pal5')
    def test_population_phase_lifetimes_and_unit_mass(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'iso.csv'
            table=pd.DataFrame(dict(Mini=[.1,.2,.4,.6,.8,.81,.82],Mass=[.1,.2,.4,.6,.8,.7,.6],int_IMF=np.arange(7.),
                Gmag=np.arange(7.),G_BPmag=np.arange(7.)+.3,G_RPmag=np.arange(7.)-.5,label=[1,1,1,1,1,3,4]))
            table.to_csv(p,index=False)
            pop=measured_population(p,0.,age=1e10,mh=-1.,source='software fixture')
            self.assertAlmostEqual(np.sum(pop.table.current_mass_msun*pop.table.number_per_msun),1.)
            # Constant MF in unevolved mass coordinate preserves tiny evolved intervals.
            self.assertLess(pop.table.number_per_msun.iloc[-1],pop.table.number_per_msun.iloc[2]/10)
            self.assertEqual(pop.metadata['weight_mode'],'empirical_present_day_ms_power_law_v1')

if __name__=='__main__':unittest.main()
