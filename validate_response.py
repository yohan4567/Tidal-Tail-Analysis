"""Model-based Poisson recovery using templates that already contain the observation response."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from fit_rate import analyse, analyse_ensemble
from provenance import sha256_file


def recovery(counts_shape,signal_draws,background,*,amplitudes=(0.,.5,1.,2.),trials=20,seed=42):
    s=np.asarray(signal_draws,float); b=np.asarray(background,float).ravel()
    if s.ndim!=2 or s.shape[1]!=b.size or trials<20 or not np.isfinite(s).all() or (s<0).any() or not np.isfinite(b).all() or (b<0).any():
        raise ValueError('Finite template ensemble, background and >=20 trials required')
    if not np.isfinite(amplitudes).all() or (np.asarray(amplitudes)<0).any():
        raise ValueError('Finite nonnegative injected amplitudes required')
    rng=np.random.default_rng(seed); rows=[]
    for amplitude in amplitudes:
        for i in range(trials):
            truth=int(rng.integers(len(s)))
            data=rng.poisson(amplitude*s[truth]+b)
            result=analyse(data,s[0],b) if len(s)==1 else analyse_ensemble(data,s,b)
            lo,hi=result['amplitude_interval_68_asymptotic']
            rows.append(dict(injected_amplitude=float(amplitude),trial=i,calibration_truth_draw=truth,
                             recovered_amplitude=result['amplitude'],interval_68_low=lo,interval_68_high=hi,
                             upper_95=result['amplitude_upper_95_one_sided_asymptotic'],covered_68=lo<=amplitude<=hi,
                             covered_upper_95=amplitude<=result['amplitude_upper_95_one_sided_asymptotic'],
                             q0_exceeds_one_sided_asymptotic_threshold=result['q0_conditional'] is not None and result['q0_conditional']>2.705543454095404))
        print(f'Recovered {trials} synthetic data sets at amplitude {amplitude:g}',flush=True)
    table=pd.DataFrame(rows)
    summary=[]
    for amplitude,t in table.groupby('injected_amplitude',sort=True):
        fields={}
        for name in ('covered_68','covered_upper_95','q0_exceeds_one_sided_asymptotic_threshold'):
            k=int(t[name].sum()); ci=binomtest(k,len(t)).proportion_ci(.95,method='exact')
            fields[name]=dict(fraction=k/len(t),interval_95=[float(ci.low),float(ci.high)])
        summary.append(dict(injected_amplitude=float(amplitude),trials=len(t),median_recovered=float(t.recovered_amplitude.median()),
                            median_bias=float(t.recovered_amplitude.median()-amplitude),**fields))
    return table,dict(experiments=summary,seed=seed,
                      interpretation='Parametric recovery under the supplied observation-response ensemble and known background shape.',
                      limitations=['Does not independently validate population/completeness/error calibration.',
                                   'Background shape and data-dependent CMD/search are not rebuilt; this is not full raw-catalogue injection.',
                                   'Asymptotic interval coverage must be assessed rather than assumed, especially at zero.',
                                   'Unknown systematics absent from the input ensemble cannot be recovered.'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--amplitudes',type=float,nargs='+',default=[0.,.5,1.,2.])
    p.add_argument('--trials',type=int,default=20)
    p.add_argument('--seed',type=int,default=42)
    a=p.parse_args()
    with np.load(a.input,allow_pickle=False) as d:
        s=d['signal_draws'] if 'signal_draws' in d else d['signal'][None,:]
        b=d['background']; n=d['counts']
    table,result=recovery(n.shape,s,b,amplitudes=a.amplitudes,trials=a.trials,seed=a.seed)
    result['input_sha256']=sha256_file(a.input)
    a.out.mkdir(parents=True,exist_ok=False)
    table.to_csv(a.out/'trials.csv',index=False)
    (a.out/'recovery.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
