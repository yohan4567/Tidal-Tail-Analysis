import unittest
import os
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).parent/'outputs/.mplcache'))

from dynamics import cylindrical_to_cartesian,cartesian_to_cylindrical,cartesian_rhs


class DynamicsChecks(unittest.TestCase):
    def test_coordinate_roundtrip(self):
        rng=np.random.default_rng(4)
        states=rng.normal(size=(50,6)); states[:,0]=rng.uniform(.2,5,50); states[:,5]=rng.uniform(-3,3,50)
        np.testing.assert_allclose(cartesian_to_cylindrical(cylindrical_to_cartesian(states)),states,atol=1e-12)

    def test_acceleration_matches_galpy_moving_potential(self):
        from galpy.orbit import Orbit
        from galpy.potential import MWPotential2014,PlummerPotential,MovingObjectPotential,evaluateRforces,evaluatezforces,evaluatephitorques
        orbit=Orbit([1.1,.03,1.0,.1,-.02,.2])
        times=np.linspace(0,-.5,501)
        orbit.integrate(times,MWPotential2014,method='dop853',rtol=1e-11,atol=1e-13)
        xyz=np.column_stack([orbit.x(times),orbit.y(times),orbit.z(times)])
        rhs=cartesian_rhs(MWPotential2014,times,xyz,.001,.01)
        pot=MWPotential2014+MovingObjectPotential(orbit=orbit,pot=PlummerPotential(amp=.001,b=.01))
        for t in (-.1013,-.25,-.4037):
            state=np.array([orbit.x(t)+.03,orbit.y(t)-.02,orbit.z(t)+.01,0.,0.,0.])
            x,y,z=state[:3]; R=np.hypot(x,y); phi=np.arctan2(y,x)
            radial=evaluateRforces(pot,R,z,phi=phi,t=t,use_physical=False)
            torque=evaluatephitorques(pot,R,z,phi=phi,t=t,use_physical=False)
            vertical=evaluatezforces(pot,R,z,phi=phi,t=t,use_physical=False)
            expected=[radial*x/R-torque*y/R**2,radial*y/R+torque*x/R**2,vertical]
            np.testing.assert_allclose(rhs(t,state)[3:],expected,rtol=2e-7,atol=1e-9)


if __name__=='__main__':
    unittest.main()
