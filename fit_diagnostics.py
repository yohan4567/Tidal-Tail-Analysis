"""Conditional count-model diagnostics with a refitted parametric bootstrap.

The background shape and selection are held fixed: this is not a calibration of
the original search, background estimation, or all astrophysical systematics.
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import xlogy,logsumexp
from scipy.stats import beta

from fit_rate import fit_ensemble,require_free_amplitude
from provenance import sha256_file


def mixture_deviance(counts,nll):
    n=np.asarray(counts,float)
    return max(0.,float(2*(np.sum(xlogy(n,n)-n)+nll)))


def diagnose(counts,signals,background,*,trials=200,seed=913):
    if int(trials)!=trials or trials<20: raise ValueError('At least 20 bootstrap trials required')
    n=np.asarray(counts,float).ravel(); s=np.asarray(signals,float); b=np.asarray(background,float).ravel()
    best=fit_ensemble(n,s,b)
    lam=best['amplitude']*s+best['background_scale']*b[None,:]
    logl=np.sum(xlogy(n[None,:],lam)-lam,axis=1)
    # Conditional mean prediction after weighting shared calibration states by
    # the entire observed dataset. Descriptive residual statistics only.
    prediction=np.exp(logl-logsumexp(logl))@lam
    residual=n-prediction
    ss=float(np.sum((n-n.mean())**2))
    descriptive=dict(ordinary_r_squared=None if ss==0 else float(1-np.sum(residual**2)/ss),
        pearson_chi_squared=float(np.sum(np.divide(residual**2,prediction,out=np.zeros_like(n),where=prediction>0))),
        fraction_expected_below_5=float(np.mean(prediction<5)),fraction_empty=float(np.mean(n==0)))
    observed=mixture_deviance(n,best['negative_log_likelihood'])
    rng=np.random.default_rng(seed); rows=[]
    for i in range(trials):
        state=int(rng.integers(len(s)))
        fake=rng.poisson(lam[state])
        fitted=fit_ensemble(fake,s,b)  # Refit BOTH amplitudes in every realization.
        rows.append(dict(trial=i,calibration_state=state,deviance=mixture_deviance(fake,fitted['negative_log_likelihood']),
            fitted_amplitude=fitted['amplitude'],fitted_background_scale=fitted['background_scale']))
        if (i+1)%20==0: print(f'Conditional goodness-of-fit bootstrap {i+1}/{trials}',flush=True)
    frame=pd.DataFrame(rows); k=int((frame.deviance>=observed).sum())
    interval=[0. if k==0 else float(beta.ppf(.025,k,trials-k+1)),
              1. if k==trials else float(beta.ppf(.975,k+1,trials-k))]
    result=dict(n_bins=len(n),n_stars=int(n.sum()),fit=best,observed_mixture_deviance=observed,
        bootstrap_trials=trials,seed=seed,exceedances=k,p_conditional_plus_one=(k+1)/(trials+1),
        p_resolution=1/(trials+1),exceedance_probability_binomial_95=interval,
        bootstrap_deviance_quantiles=dict(zip(('q025','median','q975'),map(float,np.quantile(frame.deviance,[.025,.5,.975])))),
        descriptive=descriptive,
        interpretation='Conditional parametric-bootstrap discrepancy test; small p indicates mismatch to this fixed count model, not a stream-detection significance.',
        generative_model='Draw one shared calibration state uniformly from the supplied ensemble; draw independent Poisson counts; refit both signal and background amplitudes with the same mixture likelihood.',
        limitations=['Background SHAPE was estimated from the observed field and is held fixed here, not re-estimated in each mock.',
            'Selection, CMD ridge, footprint, binning and search decisions are frozen; this is not full pipeline calibration.',
            'Only uncertainty represented in the finite calibration ensemble is generated.',
            'Ordinary R-squared and Pearson statistics are descriptive; no asymptotic chi-square p-value is reported.'])
    return result,frame,prediction


def run(input_path,out,*,trials=200,seed=913):
    input_path,out=Path(input_path),Path(out)
    if out.exists(): raise ValueError('Choose a fresh diagnostics directory')
    with np.load(input_path,allow_pickle=False) as d:
        n,b=d['counts'],d['background']
        s=d['signal_draws'] if 'signal_draws' in d else d['signal'].ravel()[None,:]
        desc=str(d['selection_description'].item())
        synthetic=bool(d['synthetic_calibration'].item()) if 'synthetic_calibration' in d else False
    try: description=json.loads(desc)
    except json.JSONDecodeError: description={}
    require_free_amplitude(description if isinstance(description,dict) else {})
    result,frame,prediction=diagnose(n,s,b,trials=trials,seed=seed)
    result.update(input_sha256=sha256_file(input_path),selection_description=desc,synthetic_calibration=synthetic)
    out.mkdir(parents=True)
    frame.to_csv(out/'bootstrap.csv',index=False)
    pd.DataFrame(dict(observed=n.ravel(),predicted=prediction)).to_csv(out/'residual_counts.csv',index=False)
    (out/'diagnostics.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parent/'outputs/.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8,4.8),layout='constrained')
    ax.hist(frame.deviance,bins=max(10,int(np.sqrt(trials))),color='#93b7df',edgecolor='white',label='Refitted model-generated counts')
    ax.axvline(result['observed_mixture_deviance'],color='#b43832',lw=2,label='Observed counts')
    ax.set(xlabel='Mixture deviance relative to saturated counts',ylabel='Bootstrap realizations',
        title=f"Conditional goodness of fit: p = {result['p_conditional_plus_one']:.3g} (B={trials})")
    ax.legend(); fig.suptitle('Fixed background shape and selection; not a full pipeline p-value',fontsize=10)
    fig.savefig(out/'goodness_of_fit.png',dpi=160); plt.close(fig)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True); p.add_argument('--out',type=Path,required=True)
    p.add_argument('--trials',type=int,default=200); p.add_argument('--seed',type=int,default=913)
    a=p.parse_args(); r=run(a.input,a.out,trials=a.trials,seed=a.seed)
    print(json.dumps({k:r[k] for k in ('observed_mixture_deviance','p_conditional_plus_one','exceedance_probability_binomial_95')},indent=2))
