"""Selected-space injection recovery conditional on a fixed background and prescribed track."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from tail_extent import measure_maps, track_vertices, project_track


def recovery(maps,track,*,inner,outer,bin_deg,half_width,true_extent,selected_stars,sigma_deg,
             trials=100,seed=42,alpha=.05,min_coverage=.7,min_consecutive=2):
    if not isinstance(trials,int) or trials<20:
        raise ValueError('At least 20 recovery trials required')
    if not np.isfinite([true_extent,selected_stars,sigma_deg]).all() or not inner<true_extent<=outer or selected_stars<=0 or sigma_deg<=0:
        raise ValueError('Require positive selected count/width and inner < true extent <= outer')
    edges=maps['edges']; bg=np.asarray(maps['background']); inside=np.asarray(maps['inside'],bool)
    c=(edges[:-1]+edges[1:])/2
    X,Y=np.meshgrid(c,c,indexing='ij')
    v,s=track_vertices(track)
    along,cross,distance=project_track(np.column_stack([X.ravel(),Y.ravel()]),v,s,(edges[1]-edges[0])/4)
    signal=np.zeros_like(along)
    for side in (-1,1):
        use=inside.ravel()&(side*along>=inner)&(side*along<true_extent)&(np.hypot(X,Y).ravel()>=inner)
        shape=np.where(use,np.exp(-.5*(cross/sigma_deg)**2),0.)
        if shape.sum()<=0:
            raise ValueError('Injected arm has no unmasked support')
        signal+=shape/shape.sum()*selected_stars/2
    signal=signal.reshape(bg.shape)
    rng=np.random.default_rng(seed)
    rows=[]
    settings=dict(inner=inner,outer=outer,bin_deg=bin_deg,half_width=half_width,alpha=alpha,
                  min_coverage=min_coverage,min_consecutive=min_consecutive)
    for i in range(trials):
        row=dict(trial=i)
        for name,mean in (('null',bg),('injected',bg+signal)):
            counts=rng.poisson(mean)
            _,result=measure_maps(counts,bg,inside,edges,track,**settings)
            row[name+'_any_candidate_segment']=any(a['candidate_segments'] for a in result.values())
            row[name+'_any_connected_extent']=any(a['connected_extent_deg'] is not None for a in result.values())
            for arm in ('backward','forward'):
                row[name+'_'+arm+'_extent_deg']=result[arm]['connected_extent_deg']
        rows.append(row)
    table=pd.DataFrame(rows)
    summary={}
    for column in ('null_any_candidate_segment','null_any_connected_extent','injected_any_candidate_segment','injected_any_connected_extent'):
        k=int(table[column].sum())
        ci=binomtest(k,trials).proportion_ci(confidence_level=.95,method='exact')
        summary[column]=dict(fraction=k/trials,n_trials=trials,interval_95=[float(ci.low),float(ci.high)])
    for arm in ('backward','forward'):
        values=table['injected_'+arm+'_extent_deg'].dropna().to_numpy()
        summary['injected_'+arm]=dict(n_recovered=len(values),
                                     median_extent_deg=float(np.median(values)) if len(values) else None,
                                     median_error_deg=float(np.median(values-true_extent)) if len(values) else None)
    summary.update(settings=settings,true_extent_deg=true_extent,selected_stars_in_unmasked_field=selected_stars,
                   sigma_deg=sigma_deg,seed=seed,
                   interpretation='Selected-space Poisson injection recovery with fixed known background.',
                   limitations=['Not an end-to-end Gaia completeness or population calibration.',
                                'Does not refit background or CMD; no spatial background mismodelling/systematics.',
                                'Single injected shape, width and amplitude; explore a grid before reporting detection limits.',
                                'Null fractions do not cover searching multiple tracks, widths, settings or clusters.'])
    return table,summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--selection',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--true-extent-deg',type=float,required=True)
    p.add_argument('--selected-stars',type=float,required=True)
    p.add_argument('--sigma-deg',type=float,default=.1)
    p.add_argument('--trials',type=int,default=100)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--bin-deg',type=float)
    p.add_argument('--half-width-deg',type=float)
    p.add_argument('--alpha',type=float,default=.05)
    p.add_argument('--min-coverage',type=float,default=.7)
    p.add_argument('--min-consecutive',type=int,default=2)
    a=p.parse_args()
    meta=json.loads((a.selection/'selection.json').read_text(encoding='utf-8'))
    with np.load(a.selection/'maps.npz',allow_pickle=False) as d:
        maps={k:d[k] for k in ('background','inside','edges')}
    pix=maps['edges'][1]-maps['edges'][0]
    table,result=recovery(maps,pd.read_csv(a.selection/'orbit_track.csv'),
                          inner=1.5*meta['tidal_radius_deg'],outer=.9*min(meta['field_radius_deg'],meta['test_radius_deg']),
                          bin_deg=a.bin_deg if a.bin_deg is not None else max(.25,3*pix),
                          half_width=a.half_width_deg if a.half_width_deg is not None else max(.25,2.1*pix),true_extent=a.true_extent_deg,
                          selected_stars=a.selected_stars,sigma_deg=a.sigma_deg,trials=a.trials,seed=a.seed,
                          alpha=a.alpha,min_coverage=a.min_coverage,min_consecutive=a.min_consecutive)
    a.out.mkdir(parents=True,exist_ok=False)
    table.to_csv(a.out/'trials.csv',index=False)
    (a.out/'recovery.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
