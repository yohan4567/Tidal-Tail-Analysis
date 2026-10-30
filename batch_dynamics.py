"""Batched Cartesian DOP853 integration with a per-tracer error controller.

Independent particles use u in [0,1], physical time t = age * (u - 1).
Batching changes evaluation order, not the force law or release distribution.
"""
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.integrate import DOP853,solve_ivp
from scipy.interpolate import CubicSpline
from dynamics import cylindrical_to_cartesian,cartesian_to_cylindrical


class TracerDOP853(DOP853):
    def _estimate_error_norm(self,K,h,scale):
        if len(scale)%6:raise ValueError('Six coordinates per tracer required')
        e5=(K.T@self.E5/scale).reshape(-1,6)
        e3=(K.T@self.E3/scale).reshape(-1,6)
        n5=np.sum(e5*e5,axis=1);n3=np.sum(e3*e3,axis=1)
        denom=np.sqrt((n5+.01*n3)*6)
        error=np.divide(abs(h)*n5,denom,out=np.zeros_like(n5),where=denom>0)
        return float(np.max(error))


def batch_rhs(mw,times,trajectory,amplitude,scale_radius,ages):
    from galpy.potential import evaluateRforces,evaluatezforces
    if any(getattr(p,'isNonAxi',False) for p in mw):raise ValueError('Axisymmetric potentials required')
    from mw2014_fast import acceleration_factory
    fast_acceleration=acceleration_factory(mw)
    times=np.asarray(times,float);order=np.argsort(times)
    spline=CubicSpline(times[order],np.asarray(trajectory)[order],axis=0,extrapolate=False)
    ages=np.asarray(ages,float);b2=float(scale_radius)**2
    if b2<=0 or not np.isfinite(ages).all() or (ages<0).any():raise ValueError('Positive scale and finite nonnegative ages required')
    if len(ages) and (-ages.max()<times.min()-1e-10 or times.max()<0):raise ValueError('Release ages outside progenitor trajectory')
    def rhs(u,flat):
        state=flat.reshape(-1,6);position=state[:,:3];r=np.hypot(position[:,0],position[:,1])
        t=ages*(u-1);centre=spline(t)
        if (r<=0).any() or not np.isfinite(centre).all():raise ValueError('Invalid time/coordinates')
        amp=np.asarray(amplitude(t) if callable(amplitude) else amplitude)
        delta=position-centre
        cluster=-np.broadcast_to(amp,(len(ages),))[:,None]*delta/(np.sum(delta*delta,axis=1)+b2)[:,None]**1.5
        force=cluster+fast_acceleration(position)
        return (np.concatenate([state[:,3:],force],axis=1)*ages[:,None]).ravel()
    return rhs


def integrate_released_batched(states,ages,mw,times,trajectory,amplitude,scale_radius,
    rtol=1e-11,atol=1e-13,checkpoint_path=None,model_identity=None,batch_size=32,checkpoint_every=100):
    initial=cylindrical_to_cartesian(states);ages=np.asarray(ages,float)
    if len(initial)!=len(ages) or not len(initial) or not np.isfinite(initial).all():raise ValueError('Finite matching nonempty states/ages required')
    if int(batch_size)!=batch_size or batch_size<1 or checkpoint_every<1:raise ValueError('Positive batch/checkpoint sizes required')
    digest=hashlib.sha256()
    for values in (initial,ages,times,trajectory,[scale_radius,rtol,atol]):
        a=np.asarray(values,float);digest.update(str(a.shape).encode());digest.update(a.tobytes())
    identity=json.dumps(dict(backend='batch_dop853_per_tracer_max_v1',model=model_identity,
                            batch_size=batch_size,inputs=digest.hexdigest()),sort_keys=True)
    checkpoint=None if checkpoint_path is None else Path(checkpoint_path);final=[]
    if checkpoint is not None and not model_identity:raise ValueError('Checkpoint requires a documented model identity')
    if checkpoint is not None and checkpoint.exists():
        with np.load(checkpoint,allow_pickle=False) as saved:
            if str(saved['identity'].item())!=identity:raise ValueError('Batch checkpoint identity mismatch')
            a=saved['final_cartesian']
            if a.ndim!=2 or a.shape[1]!=6 or len(a)>len(initial) or not np.isfinite(a).all():raise ValueError('Invalid checkpoint')
            final=list(a)
        print(f'Resuming integrated tracers: {len(final)}/{len(initial)}',flush=True)
    def save():
        if checkpoint is None:return
        temp=checkpoint.with_suffix('.tmp.npz')
        with temp.open('wb') as stream:np.savez(stream,identity=identity,final_cartesian=np.asarray(final))
        temp.replace(checkpoint)
    saved_at=len(final)
    for start in range(len(final),len(initial),batch_size):
        stop=min(start+batch_size,len(initial));a=ages[start:stop]
        rhs=batch_rhs(mw,times,trajectory,amplitude,scale_radius,a)
        solution=solve_ivp(rhs,(0.,1.),initial[start:stop].ravel(),method=TracerDOP853,rtol=rtol,atol=atol)
        if not solution.success:raise RuntimeError(solution.message)
        final.extend(solution.y[:,-1].reshape(-1,6))
        if len(final)-saved_at>=checkpoint_every:save();saved_at=len(final)
        print(f'Integrated {len(final)}/{len(initial)} tracers',flush=True)
    save()
    return cartesian_to_cylindrical(np.asarray(final))
