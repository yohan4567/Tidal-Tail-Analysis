"""Cartesian CPU integrator with a cached moving Plummer trajectory."""
import numpy as np
import hashlib
import json
from pathlib import Path
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp


def cylindrical_to_cartesian(states):
    a=np.asarray(states,float)
    R,vR,vT,z,vz,phi=np.moveaxis(a,-1,0)
    cp,sp=np.cos(phi),np.sin(phi)
    return np.stack([R*cp,R*sp,z,vR*cp-vT*sp,vR*sp+vT*cp,vz],axis=-1)


def cartesian_to_cylindrical(states):
    a=np.asarray(states,float)
    x,y,z,vx,vy,vz=np.moveaxis(a,-1,0)
    R=np.hypot(x,y)
    if (R<=0).any():
        raise ValueError('Cylindrical conversion undefined on the Galactic axis')
    return np.stack([R,(x*vx+y*vy)/R,(x*vy-y*vx)/R,z,vz,np.arctan2(y,x)],axis=-1)


def cartesian_rhs(mw,times,progenitor_xyz,amplitude,scale_radius):
    """Natural-unit acceleration: original MW forces plus analytic Plummer force.

    This path is restricted to the project's axisymmetric MW + spherical Plummer
    model. It changes numerical evaluation, not the particle release distribution.
    """
    from galpy.potential import evaluateRforces, evaluatezforces
    times=np.asarray(times,float); positions=np.asarray(progenitor_xyz,float)
    order=np.argsort(times)
    trajectory=CubicSpline(times[order],positions[order],axis=0,extrapolate=False)
    b2=float(scale_radius)**2
    if b2<=0:
        raise ValueError('Positive Plummer scale radius required')
    def rhs(t,state):
        x,y,z,vx,vy,vz=state
        R=np.hypot(x,y)
        centre=trajectory(t)
        if R<=0 or not np.isfinite(centre).all():
            raise ValueError('Orbit outside supported coordinates/time interval')
        radial=evaluateRforces(mw,R,z,t=t,use_physical=False)
        vertical=evaluatezforces(mw,R,z,t=t,use_physical=False)
        delta=state[:3]-centre
        amp=amplitude(t) if callable(amplitude) else amplitude
        cluster=-amp*delta/(np.dot(delta,delta)+b2)**1.5
        return np.array([vx,vy,vz,radial*x/R+cluster[0],radial*y/R+cluster[1],vertical+cluster[2]])
    return rhs


def integrate_released(states,ages,mw,times,trajectory,amplitude,scale_radius,rtol=1e-9,atol=1e-11,
                       checkpoint_path=None,model_identity=None,checkpoint_every=100):
    rhs=cartesian_rhs(mw,times,trajectory,amplitude,scale_radius)
    initial=cylindrical_to_cartesian(states)
    final=[]
    if len(initial)!=len(ages) or int(checkpoint_every)!=checkpoint_every or checkpoint_every<1:
        raise ValueError('Matching states/ages and a positive checkpoint interval required')
    checkpoint=None if checkpoint_path is None else Path(checkpoint_path)
    digest=hashlib.sha256()
    for array in (initial,ages,times,trajectory,np.array([scale_radius,rtol,atol])):
        a=np.asarray(array,float); digest.update(str(a.shape).encode()); digest.update(a.tobytes())
    if checkpoint is not None and not model_identity:
        raise ValueError('Checkpoints require a documented model identity')
    identity=json.dumps(dict(model=model_identity,inputs=digest.hexdigest()),sort_keys=True)
    if checkpoint is not None and checkpoint.exists():
        with np.load(checkpoint,allow_pickle=False) as saved:
            if str(saved['identity'].item())!=identity: raise ValueError('Integration checkpoint model/input mismatch')
            a=saved['final_cartesian']
            if a.ndim!=2 or a.shape[1]!=6 or len(a)>len(initial) or not np.isfinite(a).all():
                raise ValueError('Invalid checkpoint states')
            final=list(a)
        print(f'Resuming integrated tracers: {len(final)}/{len(initial)}',flush=True)
    def save_checkpoint():
        if checkpoint is None: return
        temporary=checkpoint.with_suffix('.tmp.npz')
        with temporary.open('wb') as f:
            np.savez_compressed(f,identity=identity,final_cartesian=np.asarray(final).reshape(-1,6))
        temporary.replace(checkpoint)
    for i in range(len(final),len(initial)):
        state,age=initial[i],ages[i]
        if age==0:
            final.append(state.copy())
        else:
            result=solve_ivp(rhs,(-float(age),0.),state,method='DOP853',rtol=rtol,atol=atol)
            if not result.success: raise RuntimeError(result.message)
            final.append(result.y[:,-1])
        if len(final)%checkpoint_every==0: save_checkpoint()
        if (i+1)%20==0:
            print(f'Integrated {i+1}/{len(initial)} tracers',flush=True)
    save_checkpoint()
    return cartesian_to_cylindrical(np.asarray(final))
