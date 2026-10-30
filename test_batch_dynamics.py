import os
os.environ.setdefault('MPLCONFIGDIR','outputs/.mplcache')
import unittest,tempfile
from pathlib import Path
import numpy as np
from scipy.integrate import DOP853
from batch_dynamics import batch_rhs,integrate_released_batched,TracerDOP853
from dynamics import cartesian_rhs,integrate_released,cylindrical_to_cartesian

class BatchDynamicsTests(unittest.TestCase):
    def fixture(self):
        from galpy.potential import MWPotential2014
        ts=np.linspace(-.2,0,201);xyz=np.column_stack([np.ones(201),ts,np.zeros(201)])
        states=np.array([[1.02,.01,1.,.02,.01,.01],[1.03,.02,1.,.01,.01,.02],[1.01,.01,.9,.01,0,.01]])
        return states,np.array([.1,.15,0.]),MWPotential2014,ts,xyz,lambda t:.001*(1-np.asarray(t)),.02
    def test_fast_milky_way_forces_across_scales(self):
        from galpy.potential import MWPotential2014,evaluateRforces,evaluatezforces
        from mw2014_fast import acceleration_factory
        rng=np.random.default_rng(8);xyz=rng.normal(size=(1000,3));xyz*=np.geomspace(.001,100,1000)[:,None]
        r=np.hypot(xyz[:,0],xyz[:,1]);z=xyz[:,2]
        radial=evaluateRforces(MWPotential2014,r,z,use_physical=False)
        expected=np.column_stack([radial*xyz[:,0]/r,radial*xyz[:,1]/r,evaluatezforces(MWPotential2014,r,z,use_physical=False)])
        np.testing.assert_allclose(acceleration_factory(MWPotential2014)(xyz),expected,rtol=1e-10,atol=1e-12)
        with self.assertRaises(ValueError):acceleration_factory(MWPotential2014[:2])
    def test_acceleration_and_time_rescaling(self):
        s,a,m,t,x,amp,b=self.fixture();flat=cylindrical_to_cartesian(s)
        batch=batch_rhs(m,t,x,amp,b,a)(.3,flat.ravel()).reshape(-1,6)
        scalar=cartesian_rhs(m,t,x,amp,b)
        expected=np.array([scalar(age*(.3-1),v)*age for age,v in zip(a,flat)])
        np.testing.assert_allclose(batch,expected,rtol=1e-13,atol=1e-14)
    def test_solutions_and_checkpoint(self):
        args=self.fixture();reference=integrate_released(*args)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'checkpoint.npz'
            result=integrate_released_batched(*args,batch_size=2,checkpoint_path=p,model_identity='unit')
            np.testing.assert_allclose(cylindrical_to_cartesian(result),cylindrical_to_cartesian(reference),atol=2e-8,rtol=0)
            np.testing.assert_array_equal(integrate_released_batched(*args,batch_size=2,checkpoint_path=p,model_identity='unit'),result)
            with np.load(p,allow_pickle=False) as saved:
                identity=saved['identity'].copy();prefix=saved['final_cartesian'][:2].copy()
            np.savez(p,identity=identity,final_cartesian=prefix)
            resumed=integrate_released_batched(*args,batch_size=2,checkpoint_path=p,model_identity='unit')
            np.testing.assert_allclose(resumed,result,atol=1e-13,rtol=0)
            with self.assertRaises(ValueError):integrate_released_batched(*args,batch_size=3,checkpoint_path=p,model_identity='unit')
    def test_cli_preserves_configured_backend(self):
        import sys,json
        from unittest.mock import patch
        import stream_model
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'config.json';path.write_text(json.dumps(dict(python_integrator='cartesian_batch')))
            with patch.object(sys,'argv',['stream_model.py','--config',str(path)]),patch.object(stream_model,'simulate') as call:
                stream_model.main()
                self.assertEqual(call.call_args.args[0]['python_integrator'],'cartesian_batch')
    def test_error_control_does_not_dilute_one_bad_tracer(self):
        rng=np.random.default_rng(3);K=rng.normal(size=(13,6));scale=np.ones(6)
        a=DOP853(lambda t,y:y,0,np.ones(6),1)
        single=a._estimate_error_norm(K,.1,scale)
        batch=TracerDOP853(lambda t,y:y,0,np.ones(12),1)
        together=np.column_stack([K,np.zeros_like(K)])
        self.assertAlmostEqual(batch._estimate_error_norm(together,.1,np.ones(12)),single,places=13)

if __name__=='__main__':unittest.main()
