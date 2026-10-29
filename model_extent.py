"""Weighted angular/3D radial extents of model particles, not observed detection limits."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def weighted_quantile(values,weights,q):
    values,weights=np.asarray(values,float),np.asarray(weights,float)
    if values.ndim!=1 or values.shape!=weights.shape or not len(values) or not np.isfinite(values).all() or not np.isfinite(weights).all() or (weights<0).any() or weights.sum()<=0 or not 0<q<=1:
        raise ValueError('Finite 1D values, nonnegative positive-total weights and 0 < q <= 1 required')
    good=weights>0
    values,weights=values[good],weights[good]
    order=np.argsort(values)
    index=min(np.searchsorted(np.cumsum(weights[order]),q*weights.sum(),side='left'),len(order)-1)
    return float(values[order[index]])


def summarize(particles,config,q=.95):
    p=particles
    numeric=['ra_deg','dec_deg','distance_kpc','mass_weight_msun']
    if not all(k in p for k in numeric+['arm']) or not np.isfinite(p[numeric]).all().all() or (p.distance_kpc<=0).any() or (p.mass_weight_msun<0).any():
        raise ValueError('Invalid model particle coordinates, distances or weights')
    if not set(p.arm).issubset({'leading','trailing'}):
        raise ValueError('Unrecognized model arm')
    a,d=np.radians(p.ra_deg.to_numpy()),np.radians(p.dec_deg.to_numpy())
    a0,d0=np.radians([config['ra_deg'],config['dec_deg']])
    sep=2*np.arcsin(np.sqrt(np.clip(np.sin((d-d0)/2)**2+np.cos(d)*np.cos(d0)*np.sin((a-a0)/2)**2,0,1)))
    def xyz(aa,dd,r):
        return np.array([r*np.cos(dd)*np.cos(aa),r*np.cos(dd)*np.sin(aa),r*np.sin(dd)])
    centre=xyz(a0,d0,config['distance_kpc'])
    radius3d=np.linalg.norm(xyz(a,d,p.distance_kpc.to_numpy())-centre[:,None],axis=0)*1000
    arms={}
    for name in ('leading','trailing'):
        use=p.arm.to_numpy()==name
        if not use.any():
            arms[name]=dict(status='no_particles')
            continue
        w=p.mass_weight_msun.to_numpy()[use]
        radius=weighted_quantile(np.degrees(sep[use]),w,q)
        arms[name]=dict(n_tracers=int(use.sum()),mass_msun=float(w.sum()),
                         radial_quantile_deg=radius,
                         radial_quantile_pc_at_cluster_distance=float(np.radians(radius)*config['distance_kpc']*1000),
                         radial_quantile_3d_pc=weighted_quantile(radius3d[use],w,q),
                         sampled_max_radius_deg=float(np.degrees(sep[use]).max()),
                         effective_n=float(w.sum()**2/np.sum(w*w)),
                         sampling_status='small_sample' if use.sum()<100 else 'convergence_not_established')
    return dict(cluster=config['name'],quantile=q,directions=arms,
                duration_myr=config['duration_myr'],
                interpretation='Model radial mass quantiles from the cluster centre, not along-tail length or observed extent.',
                limitations=['Depends on release duration, model, sampling and progenitor assumptions.',
                             'Model arm labels come from the release prescription.',
                             'No Gaia selection or background threshold; not comparable directly to detectable observed extent.',
                             'Physical projected scale uses D*theta at cluster distance.'])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--quantile',type=float,default=.95)
    a=p.parse_args()
    config=json.loads((a.run/'run.json').read_text(encoding='utf-8'))['config']
    result=summarize(pd.read_csv(a.run/'particles.csv'),config,a.quantile)
    a.out.mkdir(parents=True,exist_ok=False)
    (a.out/'model_extent.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
