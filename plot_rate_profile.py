"""Plot assumed mass-loss rate against conditional profiled Poisson likelihood."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from fit_rate import analyse, analyse_ensemble, ensemble_objective, profile_background,require_free_amplitude
from provenance import sha256_file

HERE=Path(__file__).resolve().parent


def scan_profile(counts,signal_draws,background,reference_rate,*,points=161,min_rate=0.,max_rate=None):
    if points<21 or int(points)!=points or not np.isfinite(reference_rate) or reference_rate<=0 or not np.isfinite(min_rate) or min_rate<0:
        raise ValueError('Need >=21 scan points, positive reference rate and nonnegative lower rate')
    s=np.asarray(signal_draws,float)
    if s.ndim!=2:
        raise ValueError('Signal draws must be [draw, bin]')
    fit=(analyse(counts,s[0],background,reference_rate) if len(s)==1 else
         analyse_ensemble(counts,s,background,reference_rate))
    best=fit['mdot_msun_myr']; lo,hi=fit['mdot_interval_68_asymptotic']; upper=fit['mdot_upper_95_one_sided_asymptotic']
    xmax=max_rate if max_rate is not None else max(2*best,1.3*upper,.01*reference_rate)
    if not np.isfinite(xmax) or xmax<=min_rate or not min_rate<=best<=xmax:
        raise ValueError('Scan range must contain the fitted best rate and have positive width')
    objective,_,n,_,b=ensemble_objective(counts,s,background)
    _,base=profile_background(objective,fit['amplitude'],n.sum()/b.sum())
    grid=np.unique(np.r_[np.linspace(min_rate,xmax,points),best,*[v for v in (lo,hi,upper) if min_rate<=v<=xmax]])
    rows=[]
    for rate in grid:
        beta,nll=profile_background(objective,rate/reference_rate,n.sum()/b.sum())
        rows.append(dict(assumed_mdot_msun_myr=float(rate),amplitude=float(rate/reference_rate),
                         fitted_background_scale=beta,profile_negative_log_likelihood=nll,
                         delta_minus2_log_likelihood=max(0.,float(2*(nll-base)))))
    return pd.DataFrame(rows),fit


def run(input_path,out,*,label=None,points=161,min_rate=0.,max_rate=None):
    input_path,out=Path(input_path),Path(out)
    with np.load(input_path,allow_pickle=False) as d:
        n,b=d['counts'],d['background']
        draws=d['signal_draws'] if 'signal_draws' in d else d['signal'].ravel()[None,:]
        reference=float(d['reference_mdot'])
        synthetic=bool(d['synthetic_calibration'].item()) if 'synthetic_calibration' in d else False
        description=str(d['selection_description'].item())
    try:
        fields=json.loads(description)
    except (json.JSONDecodeError,TypeError):
        fields={}
    if not isinstance(fields,dict): fields={}
    require_free_amplitude(fields)
    if label is None:
        for name in ('response.json','selection.json'):
            path=input_path.parent/name
            if path.exists():
                label=json.loads(path.read_text(encoding='utf-8')).get('cluster')
                if label: break
    label=label or fields.get('cluster','Stream model')
    table,fit=scan_profile(n,draws,b,reference,points=points,min_rate=min_rate,max_rate=max_rate)
    out.mkdir(parents=True,exist_ok=False)
    table.to_csv(out/'profile.csv',index=False)
    summary=dict(label=label,synthetic_calibration=synthetic,mass_basis=fields.get('mass_basis','unspecified_legacy'),
                 scientific_status=fields.get('scientific_status','conditional_model'),rate_scope=fields.get('rate_scope'),scenario=fields.get('scenario'),
                 best_mdot_msun_myr=fit['mdot_msun_myr'],interval_68_asymptotic=fit['mdot_interval_68_asymptotic'],
                 upper_95_one_sided_asymptotic=fit['mdot_upper_95_one_sided_asymptotic'],reference_mdot=reference,
                 n_calibration_draws=len(draws),input_sha256=sha256_file(input_path),
                 ordinate='q(rate) = -2 [log L_profile(rate) - log L_best]; lower is better',
                 interpretation='Relative conditional fit curve, not a p-value, posterior or absolute goodness-of-fit test.',
                 limitations=['Background normalization is refitted at every rate; background shape is held fixed.',
                              'Input shared calibration draws are marginalized jointly over all bins.',
                              'Thresholds and intervals are local/asymptotic, especially approximate near zero and for multimodal profiles.',
                              'Unknown population/calibration/model systematics and search trials are not included.'])
    (out/'profile.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    os.environ.setdefault('MPLCONFIGDIR',str(HERE/'outputs/.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    best=summary['best_mdot_msun_myr']; lo,hi=summary['interval_68_asymptotic']; upper=summary['upper_95_one_sided_asymptotic']
    x=table.assumed_mdot_msun_myr.to_numpy(); y=table.delta_minus2_log_likelihood.to_numpy()
    finite=np.isfinite(y)
    for ax in axes:
        ax.plot(x[finite],y[finite],color='#2453a6',lw=2.2,label='Profile likelihood')
        ax.axvspan(lo,hi,color='#c9ddf7',alpha=.75,label='Approx. 68% interval')
        ax.axvline(best,color='#2453a6',ls=':',lw=1.4)
        ax.scatter([best],[0],s=45,color='#2453a6',zorder=5)
        ax.axhline(1,color='#477352',ls='--',lw=1,label=r'$\Delta=1$ (~68%)')
        ax.axhline(2.705543454,color='#b06c28',ls='--',lw=1,label=r'$\Delta=2.71$ (~95% upper)')
        ax.grid(alpha=.18)
        ax.set_xlabel(r'Assumed escape rate $\dot{M}$ ($M_\odot\,\mathrm{Myr}^{-1}$)')
        ax.set_ylabel(r'$-2\,\Delta\ln\mathcal{L}$  (lower = better)')
    axes[0].set(title='Full scanned range',xlim=(x.min(),x.max()),ylim=(-.3,max(4.,float(y[finite].max())*1.05)))
    padding=max(hi-lo,.01*reference)
    axes[1].set(title=f'Best fit: {best:.4g} | 68%: [{lo:.4g}, {hi:.4g}]',
                xlim=(max(x.min(),lo-padding),min(x.max(),max(upper,hi)+.6*padding)),ylim=(-.15,7.))
    axes[1].legend(loc='upper left',fontsize=8,frameon=False)
    status='SYNTHETIC CALIBRATION — software demonstration' if synthetic else 'Conditional calibrated model'
    fig.suptitle(f'{label} | {status}',fontsize=13,color='#a62929' if synthetic else '#222222')
    fig.text(.5,.055,'Background refitted at each rate; shared calibration draws marginalized. Local/asymptotic intervals.',ha='center',fontsize=9)
    if synthetic:
        fig.text(.5,.018,'Not a measured physical mass-loss rate for this cluster.',ha='center',color='#a62929',fontsize=10)
    elif fields.get('scientific_status','').startswith('conditional_real_input'):
        fig.text(.5,.018,'Living-star mass in stated population only; mass function / extinction / zero point remain conditional.',ha='center',color='#a62929',fontsize=9)
    fig.tight_layout(rect=(0,.11,1,.93))
    fig.savefig(out/'profile.png',dpi=180)
    fig.savefig(out/'profile.svg')
    plt.close(fig)
    print(f'Profile plot: {(out/"profile.png").resolve()}',flush=True)
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--label')
    p.add_argument('--points',type=int,default=161)
    p.add_argument('--min-rate',type=float,default=0.)
    p.add_argument('--max-rate',type=float)
    a=p.parse_args()
    run(a.input,a.out,label=a.label,points=a.points,min_rate=a.min_rate,max_rate=a.max_rate)


if __name__=='__main__':
    main()
