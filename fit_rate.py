"""Conditional Poisson template fit. Requires a calibrated count template."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, minimize_scalar, brentq
from scipy.special import xlogy, logsumexp


def fit(counts, signal, background):
    """Fit nonnegative signal and background amplitudes on raw disjoint counts."""
    arrays = [np.asarray(x,dtype=float) for x in (counts,signal,background)]
    if not (arrays[0].shape == arrays[1].shape == arrays[2].shape):
        raise ValueError('Count and template shapes must match before flattening')
    n,s,b = [x.ravel() for x in arrays]
    if not (n.shape==s.shape==b.shape) or n.size<3:
        raise ValueError('Three matching arrays with at least three bins required')
    if any(not np.isfinite(x).all() or (x<0).any() for x in (n,s,b)):
        raise ValueError('Inputs must be finite and nonnegative')
    if not np.equal(n,np.round(n)).all():
        raise ValueError('Counts must be raw integer counts, not smoothed/background-subtracted')
    if s.sum()<=0 or b.sum()<=0 or np.any((s+b==0)&(n>0)):
        raise ValueError('Templates cannot explain observed counts')
    if np.linalg.matrix_rank(np.column_stack([s/s.sum(),b/b.sum()]))<2:
        raise ValueError('Signal and background shapes are not identifiable')
    def objective(v):
        lam=v[0]*s+v[1]*b
        if np.any((lam==0)&(n>0)): return np.inf
        return float(np.sum(lam-xlogy(n,lam)))
    result=minimize(objective,[max(n.sum()*.2/s.sum(),.01),max(n.sum()*.8/b.sum(),.01)],
                    bounds=[(0,None),(0,None)],method='Nelder-Mead',
                    options={'maxiter':10000,'xatol':1e-8,'fatol':1e-8})
    if not result.success: raise RuntimeError(result.message)
    return dict(amplitude=float(result.x[0]),background_scale=float(result.x[1]),
                negative_log_likelihood=float(result.fun))


def profile_summary(counts,signal,background,best=None):
    """Asymptotic profile-likelihood intervals, conditional on exact template shapes."""
    best=fit(counts,signal,background) if best is None else best
    n,s,b=[np.asarray(x,float).ravel() for x in (counts,signal,background)]
    def profile(a):
        def objective(beta):
            lam=a*s+beta*b
            if np.any((lam==0)&(n>0)):
                return np.inf
            return float(np.sum(lam-xlogy(n,lam)))
        upper=max(float(n.sum()/b.sum()),1e-12)
        r=minimize_scalar(objective,bounds=(0,upper),method='bounded',options={'xatol':1e-10})
        return min(objective(0),objective(upper),float(r.fun))
    ahat=best['amplitude']
    base=profile(ahat)
    def boundary(delta):
        f=lambda a:profile(a)-base-delta
        low=0. if f(0)<=0 else float(brentq(f,0,ahat))
        high=max(ahat+1.,1.)
        for _ in range(100):
            if f(high)>0:
                break
            high*=2
        else:
            raise RuntimeError('Could not bracket upper profile boundary')
        return [low,float(brentq(f,ahat,high))]
    null_value=profile(0)
    return dict(amplitude_interval_68_asymptotic=boundary(.5),
                amplitude_upper_95_one_sided_asymptotic=boundary(1.352771727047702)[1],
                q0_conditional=max(0.,2*(null_value-base)) if np.isfinite(null_value) else None,
                background_null_can_explain_data=bool(np.isfinite(null_value)),
                interval_caveat='Asymptotic likelihood thresholds, approximate near zero/low counts. Templates and calibration treated as exact; injection recovery required.')


def null_bootstrap(counts,signal,background,trials=200,seed=42):
    """Calibrate amplitude improvement only; does not calibrate selection/background construction."""
    if trials<20:
        raise ValueError('At least 20 null trials required')
    n,s,b=[np.asarray(x,float).ravel() for x in (counts,signal,background)]
    best=fit(n,s,b)
    def null_nll(obs):
        lam=obs.sum()/b.sum()*b
        if np.any((lam==0)&(obs>0)):
            return np.inf
        return float(np.sum(lam-xlogy(obs,lam)))
    q=max(0.,2*(null_nll(n)-best['negative_log_likelihood']))
    rng=np.random.default_rng(seed)
    mean=n.sum()/b.sum()*b
    exceed=0
    for _ in range(trials):
        draw=rng.poisson(mean)
        f=fit(draw,s,b)
        exceed+=max(0.,2*(null_nll(draw)-f['negative_log_likelihood']))>=q-1e-8
    return dict(null_trials=trials,null_seed=seed,null_p_conditional=(exceed+1)/(trials+1),
                null_p_resolution=1/(trials+1),
                null_caveat='Parametric null with fixed background/template shape and fitted null normalization; not a blind-search false-alarm probability.')


def analyse(counts,signal,background,reference_mdot=1.,null_trials=0,seed=42):
    if not np.isfinite(reference_mdot) or reference_mdot<=0:
        raise ValueError('Reference rate must be positive and finite')
    result=fit(counts,signal,background)
    result.update(profile_summary(counts,signal,background,result))
    result.update(mdot_msun_myr=result['amplitude']*reference_mdot,
                  mdot_interval_68_asymptotic=[v*reference_mdot for v in result['amplitude_interval_68_asymptotic']],
                  mdot_upper_95_one_sided_asymptotic=result['amplitude_upper_95_one_sided_asymptotic']*reference_mdot)
    if null_trials:
        result.update(null_bootstrap(counts,signal,background,null_trials,seed))
    return result


def ensemble_objective(counts,signal_draws,background):
    """Joint-data likelihood averaged over shared calibration draws (not per-bin averages)."""
    n=np.asarray(counts,float).ravel(); s=np.asarray(signal_draws,float); b=np.asarray(background,float).ravel()
    if s.ndim!=2 or s.shape[1]!=n.size or not len(s) or not np.isfinite(s).all() or (s<0).any():
        raise ValueError('signal_draws must be a nonnegative finite [calibration_draw, bin] array')
    # Reuse count/template identifiability checks and obtain a stable initialization.
    initial=fit(n,s.mean(axis=0),b)
    def objective(v):
        lam=v[0]*s+v[1]*b[None,:]
        bad=(lam==0)&(n[None,:]>0)
        logp=np.sum(xlogy(n[None,:],lam)-lam,axis=1)
        logp[np.any(bad,axis=1)]=-np.inf
        return float(-logsumexp(logp)+np.log(len(s)))
    return objective,initial,n,s,b


def fit_ensemble(counts,signal_draws,background):
    objective,initial,n,s,b=ensemble_objective(counts,signal_draws,background)
    starts=[[initial['amplitude'],initial['background_scale']],[0.,max(n.sum()/b.sum(),1e-10)]]
    answers=[minimize(objective,start,bounds=[(0,None),(0,None)],method='Nelder-Mead',
                      options={'maxiter':10000,'xatol':1e-8,'fatol':1e-8}) for start in starts]
    valid=[r for r in answers if r.success and np.isfinite(r.fun)]
    if not valid:
        raise RuntimeError('Calibration-marginal fit failed')
    best=min(valid,key=lambda r:r.fun)
    return dict(amplitude=float(best.x[0]),background_scale=float(best.x[1]),negative_log_likelihood=float(best.fun))


def profile_background(objective,amplitude,upper):
    """Return background optimum and profiled NLL for a fixed signal amplitude."""
    upper=max(float(upper),1e-12)
    grid=np.linspace(0,upper,17)
    values=np.array([objective([amplitude,v]) for v in grid])
    candidates=[(float(beta),float(value)) for beta,value in zip(grid,values)]
    for i in np.argsort(values)[:3]:
        lo,hi=grid[max(0,i-1)],grid[min(len(grid)-1,i+1)]
        if lo<hi:
            r=minimize_scalar(lambda beta:objective([amplitude,beta]),bounds=(lo,hi),method='bounded',options={'xatol':1e-10})
            candidates.append((float(r.x),float(r.fun)))
    return min(candidates,key=lambda item:item[1])


def analyse_ensemble(counts,signal_draws,background,reference_mdot=1.,null_trials=0,seed=42):
    if not np.isfinite(reference_mdot) or reference_mdot<=0:
        raise ValueError('Positive reference rate required')
    objective,initial,n,s,b=ensemble_objective(counts,signal_draws,background)
    best=fit_ensemble(n,s,b)
    def profile(a):
        return profile_background(objective,a,n.sum()/b.sum())[1]
    ahat=best['amplitude']; base=profile(ahat)
    def interval(delta):
        f=lambda a:profile(a)-base-delta
        lower=0. if f(0)<=0 else float(brentq(f,0,ahat))
        hi=max(1.,ahat+1)
        for _ in range(100):
            if f(hi)>0: break
            hi*=2
        else: raise RuntimeError('Could not bracket marginal profile interval')
        return [lower,float(brentq(f,ahat,hi))]
    i68=interval(.5); u95=interval(1.352771727047702)[1]
    best.update(amplitude_interval_68_asymptotic=i68,amplitude_upper_95_one_sided_asymptotic=u95,
                mdot_msun_myr=ahat*reference_mdot,mdot_interval_68_asymptotic=[v*reference_mdot for v in i68],
                mdot_upper_95_one_sided_asymptotic=u95*reference_mdot,n_calibration_draws=len(s),
                q0_conditional=max(0.,2*(profile(0)-base)) if np.isfinite(profile(0)) else None,
                interval_caveat='Local asymptotic profile of a finite calibration-mixture likelihood; approximate near zero/low counts and not a full confidence-set enumeration if multimodal.',
                profile_interval_scope='local_connected_region; inspect/calibrate multimodal cases separately',
                calibration_treatment='Entire-data likelihood averaged over shared input calibration draws; unrepresented systematics remain conditional.')
    if null_trials:
        if null_trials<20: raise ValueError('At least 20 null trials required')
        rng=np.random.default_rng(seed)
        mean=n.sum()/b.sum()*b
        def null_nll(obs):
            lam=obs.sum()/b.sum()*b
            return float(np.sum(lam-xlogy(obs,lam)))
        q=max(0.,2*(null_nll(n)-best['negative_log_likelihood']))
        exceed=0
        for _ in range(null_trials):
            draw=rng.poisson(mean)
            f=fit_ensemble(draw,s,b)
            exceed+=max(0.,2*(null_nll(draw)-f['negative_log_likelihood']))>=q-1e-8
        best.update(null_trials=null_trials,null_seed=seed,null_p_conditional=(exceed+1)/(null_trials+1),
                    null_p_resolution=1/(null_trials+1),null_caveat='Fixed calibration ensemble and background shape; not full search calibration.')
    return best


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,help='NPZ: counts, signal, background, reference_mdot, selection_description')
    p.add_argument('--demo',action='store_true',help='Synthetic test only; not Gaia data')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--null-trials',type=int,default=0)
    p.add_argument('--seed',type=int,default=42)
    a=p.parse_args()
    if bool(a.input)==a.demo: p.error('Choose exactly one of --input and --demo')
    if a.demo:
        x=np.linspace(-5,5,100)
        s=20*np.exp(-x*x/2); b=np.full_like(x,8.)
        n=np.random.default_rng(42).poisson(2*s+1.2*b)
        ref=1.; description='SYNTHETIC DATA: Gaussian count template, known A=2, background=1.2'
        synthetic_calibration=True
        signal_draws=None
    else:
        with np.load(a.input,allow_pickle=False) as d:
            n,s,b=[d[k] for k in ('counts','signal','background')]
            ref=float(d['reference_mdot'])
            description=str(d['selection_description'].item())
            synthetic_calibration=bool(d['synthetic_calibration'].item()) if 'synthetic_calibration' in d else False
            signal_draws=d['signal_draws'] if 'signal_draws' in d else None
        if not description.strip() or not np.isfinite(ref) or ref<=0:
            raise ValueError('Positive reference rate and documented selection required')
    try:
        description_fields=json.loads(description)
    except (json.JSONDecodeError,TypeError):
        description_fields={}
    if not isinstance(description_fields,dict):
        description_fields={}
    require_free_amplitude(description_fields)
    result=(analyse(n,s,b,ref,a.null_trials,a.seed) if signal_draws is None else
            analyse_ensemble(n,signal_draws,b,ref,a.null_trials,a.seed))
    result.update(mdot_msun_myr=result['amplitude']*ref,reference_mdot=ref,
                  selection_description=description,synthetic=a.demo,
                  synthetic_calibration=synthetic_calibration,
                  mass_basis=description_fields.get('mass_basis','unspecified_legacy'),
                  rate_scope=description_fields.get('rate_scope','Mass represented by supplied population only'),
                  scientific_status=description_fields.get('scientific_status','conditional_model'),scenario=description_fields.get('scenario'),
                  interpretation='Conditional rate and asymptotic statistical intervals. No unconditional detection claim.',
                  limitation='Requires expected SELECTED STAR COUNTS at reference rate; raw tracer mass is NOT a count template.')
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('x',encoding='utf-8') as f: json.dump(result,f,indent=2)
    print(json.dumps(result,indent=2))


def require_free_amplitude(description):
    if description.get('amplitude_policy')=='fixed_reference_rate' or description.get('dynamics_mass_history') is not None:
        raise ValueError('This template has rate-coupled dynamics. Use score_coupled_grid.py; changing a free amplitude is not a physical rate fit.')


if __name__=='__main__': main()
