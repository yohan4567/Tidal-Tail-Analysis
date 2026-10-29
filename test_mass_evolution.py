import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).parent/'outputs/.mplcache'))
import numpy as np
import pandas as pd

from mass_evolution import MassEvolution,from_config
from dynamics import cartesian_rhs,integrate_released
from fit_rate import require_free_amplitude
from response_model import history_weights
from score_coupled_grid import score


class EvolvingMassChecks(unittest.TestCase):
    def test_mass_budget_and_domain(self):
        e=MassEvolution(10000,1000,5,2,'explicit test budget')
        np.testing.assert_allclose(e.mass([0,500,1000]),[10000,15000,20000])
        self.assertAlmostEqual(e.ratio(1000),2.)
        with self.assertRaises(ValueError): e.mass(1001)
        with self.assertRaises(ValueError): MassEvolution(10000,1000,5,-1,'invalid')
        np.testing.assert_allclose(MassEvolution(10000,1000,0,1,'constant limit').ratio([0,1000]),1)

    def test_explicit_mass_basis_and_free_amplitude_guard(self):
        self.assertIsNone(from_config({}))
        c=dict(mass_msun=10000,duration_myr=1000,reference_mdot_msun_myr=5,
            mass_evolution=dict(mode='linear_escape',dynamical_to_population_mass=1,source_description='test'))
        with self.assertRaises(ValueError): from_config(c)
        c['mass_evolution']['population_mass_basis']='isochrone_current_stellar_mass'
        self.assertAlmostEqual(from_config(c).mass(1000),15000)
        with self.assertRaisesRegex(ValueError,'rate-coupled'): require_free_amplitude(dict(amplitude_policy='fixed_reference_rate'))
        with self.assertRaises(ValueError): history_weights(pd.DataFrame(),dict(dynamics_mass_history={'mode':'test'}),pd.DataFrame())

    def test_release_distribution_uses_instantaneous_mass(self):
        from galpy.df import chen24spraydf
        from galpy.orbit import Orbit
        from galpy.potential import MWPotential2014,PlummerPotential
        original=Orbit.integrate
        def integrate(self,t,pot,*a,**kw):
            kw.update(method='dop853',rtol=1e-10,atol=1e-12,progressbar=False)
            return original(self,t,pot,*a,**kw)
        def make(mass):
            return chen24spraydf(mass,progenitor=Orbit([1.1,0.,1.,.1,0.,0.]),pot=MWPotential2014,
                progpot=PlummerPotential(amp=1e-5,b=.001),tdisrupt=.5,tail='both')
        with patch.object(Orbit,'integrate',integrate):
            fixed=make(1e-5); constant=make(lambda t: 1e-5*np.ones_like(t)); evolving=make(lambda t:1e-5*(1-2*np.asarray(t)))
        xyz=np.array([[1.1,0,.1]]*3); vel=np.array([[0,1.,0]]*3); dt=np.array([.1,.2,.3])
        def release(df):
            np.random.seed(7); return np.array(df.spray_df(xyz,vel,dt)).T
        a,b,c=release(fixed),release(constant),release(evolving)
        np.testing.assert_array_equal(a,b)
        baseline=np.column_stack([xyz,vel]); ratio=(1+2*dt)**(1/3)
        np.testing.assert_allclose(c-baseline,(a-baseline)*ratio[:,None],rtol=2e-10,atol=1e-12)

    def test_evolving_force_matches_galpy_wrapper(self):
        from galpy.orbit import Orbit
        from galpy.potential import MWPotential2014,PlummerPotential,TimeDependentAmplitudeWrapperPotential,MovingObjectPotential,evaluateRforces,evaluatezforces,evaluatephitorques
        orbit=Orbit([1.1,.03,1.,.1,-.02,.2]); times=np.linspace(0,-.5,501)
        orbit.integrate(times,MWPotential2014,method='dop853',rtol=1e-11,atol=1e-13)
        xyz=np.column_stack([orbit.x(times),orbit.y(times),orbit.z(times)])
        rhs=cartesian_rhs(MWPotential2014,times,xyz,lambda t:.001*(1-2*t),.01)
        moving=MovingObjectPotential(orbit=orbit,pot=TimeDependentAmplitudeWrapperPotential(pot=PlummerPotential(amp=.001,b=.01),A=lambda t:1-2*t))
        pot=MWPotential2014+moving
        for t in (-.1013,-.25,-.4037):
            state=np.array([orbit.x(t)+.03,orbit.y(t)-.02,orbit.z(t)+.01,0.,0.,0.]); x,y,z=state[:3]
            R=np.hypot(x,y); phi=np.arctan2(y,x)
            rf=evaluateRforces(pot,R,z,phi=phi,t=t,use_physical=False); torque=evaluatephitorques(pot,R,z,phi=phi,t=t,use_physical=False)
            expected=[rf*x/R-torque*y/R**2,rf*y/R+torque*x/R**2,evaluatezforces(pot,R,z,phi=phi,t=t,use_physical=False)]
            np.testing.assert_allclose(rhs(t,state)[3:],expected,rtol=2e-7,atol=1e-9)

    def test_checkpoint_resumes_without_reintegrating_and_rejects_changes(self):
        states=np.array([[1.,0.,1.,.1,0.,0.],[1.1,0.,1.,.1,0.,.1]])
        ages=np.array([.1,.2]); times=np.array([-1.,0.]); xyz=np.zeros((2,3))
        done=SimpleNamespace(success=True,y=np.array([[1.],[0.],[.1],[0.],[1.],[0.]]))
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent/'outputs') as tmp:
            f=Path(tmp)/'checkpoint.npz'
            args=(states,ages,[],times,xyz,.001,.01)
            with patch('dynamics.solve_ivp',side_effect=[done,RuntimeError('interrupted')]):
                with self.assertRaises(RuntimeError): integrate_released(*args,checkpoint_path=f,model_identity='test',checkpoint_every=1)
            with patch('dynamics.solve_ivp',return_value=done) as solver:
                result=integrate_released(*args,checkpoint_path=f,model_identity='test',checkpoint_every=1)
                self.assertEqual(solver.call_count,1); self.assertEqual(result.shape,(2,6))
            with self.assertRaises(ValueError): integrate_released(*args,checkpoint_path=f,model_identity='different')

    def test_coupled_grid_fits_background_not_signal_amplitude(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); paths=[]
            for rate in (1.,2.):
                p=root/f'{rate}.npz'; paths.append(p)
                desc=dict(observation_identity='same',selection_sha256='same',mass_basis='isochrone_current_stellar_mass',
                    amplitude_policy='fixed_reference_rate',dynamics_mass_history=dict(stellar_rate_msun_myr=rate),scenario=dict(escape_mass_tilt=0))
                np.savez(p,counts=[15,7,5],background=[1.,1.,1.],signal=np.array([5.,1.,0.])*rate,
                    reference_mdot=rate,selection_description=json.dumps(desc))
            result=score(paths,root/'scores')
            rows=result['candidates']
            self.assertEqual([r['signal_amplitude'] for r in rows],[1.,1.])
            self.assertLess(rows[1]['negative_log_likelihood'],rows[0]['negative_log_likelihood'])


if __name__=='__main__': unittest.main()
