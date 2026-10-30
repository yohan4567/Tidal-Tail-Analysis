import unittest,tempfile,json,hashlib
from pathlib import Path
from importlib.metadata import version
import numpy as np
from grid_reuse import reuse_simulation
from provenance import sha256_file

class GridReuseTests(unittest.TestCase):
    def fixture(self,root,complete):
        src=root/'source';src.mkdir();config=dict(n_particles=2,seed=5)
        signature=dict(config=config,source_sha256={'batch_dynamics.py':sha256_file('batch_dynamics.py')},versions={'numpy':version('numpy')})
        (src/'simulation_input.json').write_text(json.dumps(signature))
        if complete:
            (src/'particles.csv').write_text('value\n1\n2\n')
            (src/'run.json').write_text(json.dumps(dict(config=config,particles_sha256=sha256_file(src/'particles.csv'))))
        else:
            model=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()
            np.savez(src/'integration_checkpoint.npz',identity=json.dumps(dict(model=model,inputs='validated_by_integrator_on_resume')),final_cartesian=np.ones((1,6)))
        return src,config
    def test_completed_output_reused_and_wrong_config_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            r=Path(tmp);src,c=self.fixture(r,True)
            a=reuse_simulation(src,r/'target/simulation',c)
            self.assertTrue(a['complete']);self.assertEqual(a['retained_tracers'],2)
            with self.assertRaises(ValueError):reuse_simulation(src,r/'bad',dict(c,seed=6))
    def test_partial_checkpoint_and_corruption_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            r=Path(tmp);src,c=self.fixture(r,False)
            a=reuse_simulation(src,r/'target/simulation',c)
            self.assertFalse(a['complete']);self.assertEqual(a['retained_tracers'],1)
            np.savez(src/'integration_checkpoint.npz',identity=json.dumps(dict(model='wrong')),final_cartesian=np.ones((1,6)))
            with self.assertRaises(ValueError):reuse_simulation(src,r/'bad',c)
    def test_missing_source_is_new_work(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertIsNone(reuse_simulation(Path(tmp)/'absent',Path(tmp)/'new',{}))

if __name__=='__main__':unittest.main()
