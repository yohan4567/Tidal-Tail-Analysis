"""Check resource controls independently of numerical and observation settings."""
import json,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import portable_cluster as portable


class WorkerTests(unittest.TestCase):
    def test_dynamics_and_response_limits_are_independent(self):
        with tempfile.TemporaryDirectory() as d:
            base=Path(d);out=base/'results';out.mkdir()
            (base/'portable_manifest.json').write_text('{}')
            for offset in ('offset_0','offset_minus0p5','offset_minus1'):
                folder=out/offset;folder.mkdir();(folder/'population_prior.json').write_text(json.dumps({'response_configs':['example.json']}))
            calls=[]
            def grid(*args,**kwargs):calls.append(kwargs)
            modules={'run_joint_study':SimpleNamespace(run=grid),'bayesian_grid':SimpleNamespace(run=lambda *args:None),
                'production_study':SimpleNamespace(RATES=[1.,3.],TIMES=[100.,200.],SEEDS=[1009])}
            with patch.dict('sys.modules',modules):
                portable.run(base,out,{'cluster':'NGC362','cluster_age_myr':10750.},24,8)
            self.assertEqual(len(calls),3)
            self.assertTrue(all(c['workers']==24 and c['response_workers']==8 for c in calls))
            self.assertTrue(all(c['rates']==[1.,3.] and c['times']==[100.,200.] for c in calls))
            settings=json.loads((out/'execution_settings.json').read_text())
            self.assertEqual(settings['workers'],24);self.assertEqual(settings['response_workers'],8)

    def test_existing_package_can_be_used_by_updated_clone(self):
        with patch('sys.argv',['portable_cluster.py','status','--bundle','example','--workers','24']),patch.object(portable,'status') as status:
            portable.main()
            self.assertEqual(status.call_args.args[0],Path('example').resolve()/'results')

    def test_nonpositive_counts_rejected(self):
        for args in (['--workers','0'],['--response-workers','0']):
            with patch('sys.argv',['portable_cluster.py','status',*args]),self.assertRaises(SystemExit):portable.main()


if __name__=='__main__':unittest.main()
