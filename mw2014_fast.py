"""Algebraically equivalent vector forces for the three MWPotential2014 components.

Parameters are read from galpy objects; unsupported potential classes are refused.
"""
import numpy as np
from scipy.special import gamma,gammainc


def acceleration_factory(potentials):
    from galpy.potential import PowerSphericalPotentialwCutoff,MiyamotoNagaiPotential,NFWPotential
    classes=(PowerSphericalPotentialwCutoff,MiyamotoNagaiPotential,NFWPotential)
    if len(potentials)!=3 or {type(p) for p in potentials}!=set(classes):raise ValueError('Exactly the MWPotential2014 component classes required')
    bulge=next(p for p in potentials if type(p) is classes[0]);disk=next(p for p in potentials if type(p) is classes[1]);halo=next(p for p in potentials if type(p) is classes[2])
    k=(3-bulge.alpha)/2
    if k<=0 or bulge.rc<=0 or disk._b<=0 or halo.a<=0:raise ValueError('Unsupported potential parameters')
    bulge_scale=bulge._amp*2*np.pi*bulge.rc**(3-bulge.alpha)*gamma(k)
    rc2=bulge.rc**2;da=disk._a;db2=disk._b**2;disk_amp=disk._amp;ha=halo.a;halo_amp=halo._amp
    def acceleration(xyz):
        xyz=np.asarray(xyz,float);r2=np.sum(xyz*xyz,axis=-1);r=np.sqrt(r2)
        if (r<=0).any():raise ValueError('Galactic origin is unsupported')
        mass_bulge=bulge_scale*gammainc(k,r2/rc2)
        spherical=(-mass_bulge+halo_amp*(r/(ha+r)-np.log1p(r/ha)))/(r2*r)
        result=spherical[...,None]*xyz
        z=xyz[...,2];zb=np.sqrt(z*z+db2);az=da+zb
        disk_force=-disk_amp/(xyz[...,0]**2+xyz[...,1]**2+az*az)**1.5
        result[...,:2]+=disk_force[...,None]*xyz[...,:2]
        result[...,2]+=disk_force*z*az/zb
        return result
    return acceleration
