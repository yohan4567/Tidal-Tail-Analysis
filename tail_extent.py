"""Operational projected extent along a predeclared track, conditional on a background map."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import poisson

HERE=Path(__file__).resolve().parent


def track_vertices(track):
    """Orient a two-sided track from negative-side end through centre to positive end."""
    arrays=[]
    for side in (-1,0,1):
        v=track.loc[track.side==side,['x_deg','y_deg']].to_numpy(float)
        if side==-1:
            v=v[::-1]
        arrays.append(v)
    if len(arrays[1])!=1 or not len(arrays[0]) or not len(arrays[2]):
        raise ValueError('Track requires one centre (side=0), and outward-ordered side=-1/+1 points')
    v=np.concatenate(arrays)
    if not np.isfinite(v).all():
        raise ValueError('Nonfinite track vertices')
    ds=np.linalg.norm(np.diff(v,axis=0),axis=1)
    s=np.r_[0.,np.cumsum(ds)]
    s-=s[len(arrays[0])]
    # Consecutive duplicates can occur at t=0; remove them without shifting the origin.
    keep=np.r_[True,ds>1e-10]
    return v[keep],s[keep]


def project_track(points,vertices,s,spacing):
    """Densify, then project to either segment adjacent to the nearest vertex."""
    if spacing<=0 or len(vertices)<3:
        raise ValueError('Positive spacing and at least 3 track vertices required')
    dense_s=np.unique(np.r_[np.arange(s[0],s[-1],spacing),s,0.])
    v=np.column_stack([np.interp(dense_s,s,vertices[:,i]) for i in range(2)])
    nearest=cKDTree(v).query(points)[1]
    best=np.full(len(points),np.inf)
    along=np.zeros(len(points)); cross=np.zeros(len(points))
    for shift in (-1,0):
        k=np.clip(nearest+shift,0,len(v)-2)
        direction=v[k+1]-v[k]
        norm2=np.einsum('ij,ij->i',direction,direction)
        t=np.clip(np.einsum('ij,ij->i',points-v[k],direction)/norm2,0,1)
        delta=points-(v[k]+t[:,None]*direction)
        dist=np.linalg.norm(delta,axis=1)
        better=dist<best
        best[better]=dist[better]
        along[better]=(dense_s[k]+t*(dense_s[k+1]-dense_s[k]))[better]
        signed=(direction[:,0]*delta[:,1]-direction[:,1]*delta[:,0])/np.sqrt(norm2)
        cross[better]=signed[better]
    return along,cross,best


def measure_maps(counts,background,inside,edges,track,*,inner,outer,bin_deg,half_width,
                 alpha=.05,min_coverage=.7,min_consecutive=2):
    """Disjoint pixel bins, Bonferroni local Poisson tests over both track directions.

    The supplied background is treated as known. Results are conditional candidates,
    not a calibrated all-sky detection probability or a physical stream endpoint.
    """
    counts,background=np.asarray(counts,float),np.asarray(background,float)
    inside=np.asarray(inside,bool)
    if counts.shape!=background.shape or counts.shape!=inside.shape or counts.shape!=(len(edges)-1,len(edges)-1):
        raise ValueError('Map shapes and grid edges do not match')
    if not np.isfinite(counts).all() or not np.isfinite(background).all() or (counts<0).any() or (background<0).any() or not np.equal(counts,np.round(counts)).all():
        raise ValueError('Raw nonnegative integer counts and nonnegative finite background required')
    steps=np.diff(edges)
    if not np.isfinite(edges).all() or not (steps>0).all() or not np.allclose(steps,steps[0]):
        raise ValueError('Finite uniformly spaced increasing grid required')
    pix=steps[0]
    if not np.isfinite([inner,outer,bin_deg,half_width,alpha,min_coverage]).all() or not (0<=inner<outer) or bin_deg<pix or half_width<2*pix:
        raise ValueError('Require 0 <= inner < outer, bin >= 1 pixel and half-width >= 2 pixels')
    if not 0<alpha<1 or not 0<min_coverage<=1 or min_consecutive<1 or int(min_consecutive)!=min_consecutive:
        raise ValueError('Invalid alpha, coverage or consecutive-bin requirement')
    bounds=np.arange(inner,outer,bin_deg)
    if len(bounds)<2:
        raise ValueError('Tail region must contain at least two complete bins')
    # Omit a very short terminal bin; retain one at least half the requested width.
    bounds=np.r_[bounds,outer]
    if bounds[-1]-bounds[-2]<.5*bin_deg:
        bounds=bounds[:-1]
    v,s=track_vertices(track)
    cen=(edges[:-1]+edges[1:])/2
    X,Y=np.meshgrid(cen,cen,indexing='ij')
    points=np.column_stack([X.ravel(),Y.ravel()])
    along,cross,distance=project_track(points,v,s,pix/4)
    radial=np.hypot(points[:,0],points[:,1])
    n,b,ok=counts.ravel(),background.ravel(),inside.ravel()
    rows=[]
    for side,name in ((-1,'backward'),(1,'forward')):
        for lo,hi in zip(bounds[:-1],bounds[1:]):
            # All counts stay in a single pixel/bin; the map is never smoothed.
            geometry=(side*along>=lo)&(side*along<hi)&(distance<=half_width)&(radial>=inner)
            use=geometry&ok
            area=use.sum()*pix*pix
            # The full strip is needed to flag a track leaving the downloaded grid.
            ideal=2*half_width*(hi-lo)
            coverage=min(area/max(ideal,1e-30),1.)
            in_track=hi<=max(side*s)
            eligible=coverage>=min_coverage and use.sum()>=4 and in_track
            observed=int(n[use].sum()); expected=float(b[use].sum())
            excess=observed-expected
            second=float(np.sum((n[use]-b[use])*cross[use]**2))
            first=float(np.sum((n[use]-b[use])*cross[use]))
            var=second/excess-(first/excess)**2 if excess>0 else -1
            width=float(np.sqrt(var)) if var>=0 else None
            rows.append(dict(direction=name,lo_deg=float(lo),hi_deg=float(hi),mid_deg=float((lo+hi)/2),
                             observed=observed,background=expected,excess=excess,
                             area_deg2=float(area),coverage=float(coverage),n_pixels=int(use.sum()),
                             eligible=bool(eligible),p_local=float(poisson.sf(observed-1,expected)) if eligible else None,
                             aperture_rms_width_deg=width,transverse_excess_mean_deg=float(first/excess) if excess>0 else None))
    tests=sum(r['eligible'] for r in rows)
    for r in rows:
        r['p_bonferroni']=min(1.,r['p_local']*tests) if r['eligible'] else None
        r['candidate']=bool(r['eligible'] and r['excess']>0 and r['p_bonferroni']<=alpha)
        if not r['candidate'] or r['aperture_rms_width_deg'] is None or r['aperture_rms_width_deg']>half_width:
            r['aperture_rms_width_deg']=None
    summary={}
    for name in ('backward','forward'):
        arm=[r for r in rows if r['direction']==name]
        runs=[]; current=[]
        for i,r in enumerate(arm):
            if r['candidate']:
                current.append(i)
            elif current:
                if len(current)>=min_consecutive: runs.append(current)
                current=[]
        if len(current)>=min_consecutive: runs.append(current)
        segments=[dict(start_deg=arm[g[0]]['lo_deg'],end_deg=arm[g[-1]]['hi_deg'],n_bins=len(g)) for g in runs]
        connected=next((g for g in runs if g[0]==0),None)
        limit=None
        if connected:
            end=connected[-1]+1
            limit=('field_or_search_boundary' if end==len(arm) else
                   'coverage_or_track_limit' if not arm[end]['eligible'] else 'below_conditional_threshold')
        summary[name]=dict(candidate_segments=segments,
                           connected_extent_deg=arm[connected[-1]]['hi_deg'] if connected else None,
                           extent_bracket_deg=[arm[connected[-1]]['lo_deg'],arm[connected[-1]]['hi_deg']] if connected else None,
                           limiting_reason=limit,status='conditional_candidate' if connected else 'no_connected_candidate',
                           note='No connected candidate is not evidence that no tail exists.')
    return pd.DataFrame(rows),summary


def run(selection,out,bin_deg=None,half_width=None,alpha=.05,min_coverage=.7,min_consecutive=2):
    selection,out=Path(selection),Path(out)
    report=json.loads((selection/'selection.json').read_text(encoding='utf-8'))
    required=('tidal_radius_deg','distance_kpc','field_radius_deg')
    if any(k not in report for k in required) or not (selection/'orbit_track.csv').exists():
        raise ValueError('Re-run selection_bridge.py to export geometry/track metadata')
    track=pd.read_csv(selection/'orbit_track.csv')
    with np.load(selection/'maps.npz',allow_pickle=False) as d:
        maps={k:d[k] for k in ('counts','background','inside','edges')}
    pix=float(maps['edges'][1]-maps['edges'][0])
    bin_deg=bin_deg if bin_deg is not None else max(.25,3*pix)
    half_width=half_width if half_width is not None else max(.25,2.1*pix)
    inner=1.5*report['tidal_radius_deg']
    outer=.9*min(report['field_radius_deg'],report.get('test_radius_deg',5.))
    table,arms=measure_maps(**maps,track=track,inner=inner,outer=outer,bin_deg=bin_deg,
                            half_width=half_width,alpha=alpha,min_coverage=min_coverage,min_consecutive=min_consecutive)
    factor=report['distance_kpc']*1000*np.pi/180
    for arm in arms.values():
        length=arm['connected_extent_deg']
        arm['projected_extent_pc_at_cluster_distance']=length*factor if length is not None else None
    table['aperture_rms_width_pc']=table.aperture_rms_width_deg*factor
    result=dict(cluster=report['cluster'],directions=arms,bin_deg=bin_deg,half_width_deg=half_width,
                inner_exclusion_deg=inner,outer_search_deg=outer,alpha=alpha,min_coverage=min_coverage,
                min_consecutive=min_consecutive,n_tested_bins=int(table.eligible.sum()),
                interpretation='Conditional orbit-aligned candidate extent; not a physical tail radius or a globally calibrated detection.',
                limitations=['Same-field fitted background treated as known; shape/normalization uncertainty not propagated.',
                             'Orbit-following membership cuts can bias orbit-aligned structure.',
                             'Bonferroni covers only these two-direction bins, not scans over width, direction, settings or clusters.',
                             'Forward/backward orbit directions are not a confirmed leading/trailing identification.',
                             'Lengths use tangent-plane arc at the cluster distance; no 3D deprojection.',
                             'Widths are background-subtracted second moments within the chosen aperture, not deconvolved Gaussian widths.',
                             'Count error bars are exploratory sqrt(N+B); they are not confidence intervals for fitted-background uncertainty.',
                             'Extent brackets describe bin resolution, not confidence intervals; injection recovery is required.'])
    out.mkdir(parents=True,exist_ok=False)
    table.to_csv(out/'profile.csv',index=False)
    (out/'extent.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    os.environ.setdefault('MPLCONFIGDIR',str(HERE/'outputs/.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,3,figsize=(15,4.5),constrained_layout=True)
    edges=maps['edges']; H=maps['counts']; inside=maps['inside']
    image=axs[0].imshow(np.where(inside,H,np.nan).T,origin='lower',extent=[edges[0],edges[-1],edges[0],edges[-1]],cmap='viridis')
    for side,color,name in ((-1,'orange','Backward orbit'),(1,'cyan','Forward orbit')):
        armtrack=track[(track.side==0)|(track.side==side)]
        axs[0].plot(armtrack.x_deg,armtrack.y_deg,color=color,lw=1,label=name)
    axs[0].set(xlabel='Gnomonic east (deg)',ylabel='Gnomonic north (deg)',title='Observed candidates and prescribed orbit',xlim=(edges[0],edges[-1]),ylim=(edges[0],edges[-1]))
    axs[0].legend(fontsize=7); fig.colorbar(image,ax=axs[0],label='Stars / pixel')
    for name,color in (('backward','tab:orange'),('forward','tab:blue')):
        t=table[table.direction==name]
        axs[1].errorbar(t.mid_deg,t.excess,yerr=np.sqrt(t.observed+t.background),marker='.',color=color,label=name)
        detected=t[t.candidate]
        axs[1].scatter(detected.mid_deg,detected.excess,facecolors='none',edgecolors='red',s=80)
        axs[2].plot(t.mid_deg,t.aperture_rms_width_deg,marker='.',color=color,label=name)
    axs[1].axhline(0,color='gray',ls='--'); axs[1].legend()
    axs[1].set(xlabel='Along-track distance (tangent deg)',ylabel='Observed - background counts',title='Red circles: conditional candidate bins')
    axs[2].set(xlabel='Along-track distance (tangent deg)',ylabel='Aperture RMS width (deg)',title='Exploratory width; no PSF/error correction')
    if not table.aperture_rms_width_deg.notna().any():
        axs[2].text(.5,.5,'No resolved candidate widths',ha='center',va='center',transform=axs[2].transAxes)
    fig.suptitle(f"{report['cluster']} | prescribed orbit | conditional on fitted background")
    fig.savefig(out/'extent.png',dpi=160); plt.close(fig)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--selection',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--bin-deg',type=float)
    p.add_argument('--half-width-deg',type=float)
    p.add_argument('--alpha',type=float,default=.05)
    p.add_argument('--min-coverage',type=float,default=.7)
    p.add_argument('--min-consecutive',type=int,default=2)
    a=p.parse_args()
    result=run(a.selection,a.out,a.bin_deg,a.half_width_deg,a.alpha,a.min_coverage,a.min_consecutive)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
