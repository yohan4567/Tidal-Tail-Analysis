"""Inspect a physical population against core stars BEFORE applying the CMD cut."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from population_model import Population,build_population,read_isochrone
from provenance import sha256_file
from selection_bridge import DEFAULT_STUDY,load_study


def core_sample(selection,radius=.15):
    selection=Path(selection)
    report=json.loads((selection/'selection.json').read_text(encoding='utf-8'))
    T=load_study(DEFAULT_STUDY,report['legacy'])
    row=report['cluster_selection_parameters']; field=Path(report['source_field'])
    if sha256_file(field)!=report['source_field_sha256']: raise ValueError('Observation changed')
    raw=pd.read_csv(field,low_memory=False)
    d=T.quality(raw).copy()
    x,y=T.gnomonic(d.ra.to_numpy(),d.dec.to_numpy(),row['ra'],row['dec'])
    keep=(np.hypot(x,y)<radius)&T.plx_select(d,row['R_sun'])&T.pm_select(d,(row['pmra'],row['pmdec']),row['R_sun'],report['sig0_kms'])
    d=d.loc[keep]
    return d,report,field,T


def run(population,selection,out,*,radius=.15,ebv=.03):
    selection,out=Path(selection),Path(out)
    pop=Population.load(population)
    d,report,field,T=core_sample(selection,radius)
    row=report['cluster_selection_parameters']
    av=3.1*ebv; ag,abp,arp=np.array([.83627,1.08337,.63439])*av
    t=pop.table; colour=t.M_BP-t.M_RP+abp-arp
    g=t.M_G+5*np.log10(row['R_sun'])+10+ag
    os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parent/'outputs/.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,7),layout='constrained')
    ax.scatter(d.bp_rp,d.phot_g_mean_mag,s=9,c='#777777',alpha=.45,label=f'Core data before CMD cut (N={len(d)})')
    phase=pd.to_numeric(t.phase,errors='coerce')
    for label,mask in [('Main sequence / giant branch',phase<=3),('Later evolutionary phases',(phase>3)&(phase<9))]:
        ax.scatter(colour[mask],g[mask],s=11,label=f'PARSEC: {label}')
    ax.set(xlabel='Gaia BP - RP',ylabel='Gaia G (mag)',xlim=(-.6,2.5),ylim=(21,11.5),
        title=f"{report['cluster']}: physical population / observed core CMD")
    ax.legend(fontsize=8); ax.grid(alpha=.2)
    fig.suptitle('Adopted age/abundance/distance/extinction; inspect agreement, not a fitted age',fontsize=10)
    out.mkdir(parents=True,exist_ok=False)
    fig.savefig(out/'population_cmd.png',dpi=160); plt.close(fig)
    meta=dict(n_core_stars=len(d),radius_deg=radius,cmd_cut_applied=False,population=pop.metadata,
        source_field=str(field),source_sha256=sha256_file(field),ebv=ebv,distance_kpc=row['R_sun'],
        interpretation='Visual conditional CMD check, not an independent fit or validation of the escape mass function. Field contamination and binaries remain.')
    (out/'cmd_check.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    return meta


def compare_isochrones(isochrones,selection,out,*,radius=.15,ebv=.03):
    """Diagnostic HB sensitivity: no automatic membership or eta/age inference."""
    out=Path(out)
    if out.exists(): raise ValueError('Use a fresh HB comparison directory')
    d,report,field,T=core_sample(selection,radius)
    row=report['cluster_selection_parameters']; distance=row['R_sun']
    av=3.1*ebv; ag,abp,arp=np.array([.83627,1.08337,.63439])*av
    box=dict(g_min=16.8,g_max=17.8,color_min=-.5,color_max=.85)
    def in_box(g,c):
        return (g>=box['g_min'])&(g<=box['g_max'])&(c>=box['color_min'])&(c<=box['color_max'])
    observed_box=int(in_box(d.phot_g_mean_mag,d.bp_rp).sum())
    ridge=report['cmd_ridge']
    observed_cmd=T.cmd_select(d,np.array(ridge['G']),np.array(ridge['color']),np.array(ridge['width']))
    observed_box_retained=int((in_box(d.phot_g_mean_mag,d.bp_rp)&observed_cmd).sum())
    models=[]
    for path in map(Path,isochrones):
        imeta=json.loads(path.with_suffix('.json').read_text(encoding='utf-8'))
        if sha256_file(path)!=imeta['sha256'] or imeta.get('av')!=0 or imeta.get('synthetic_calibration'):
            raise ValueError('HB comparison needs verified unextinguished physical isochrones')
        pop=build_population(read_isochrone(path),mass_basis='isochrone_current_stellar_mass',source_description=imeta['source_description'],
            age_years=imeta['age_years'],metallicity_dex=imeta['metallicity_mh'])
        eta=float(imeta.get('eta_reimers',imeta['request_parameters']['eta_reimers']))
        t=pop.table; colour=t.M_BP-t.M_RP+abp-arp; g=t.M_G+5*np.log10(distance)+10+ag
        phase=pd.to_numeric(t.phase); hb=(phase==4).to_numpy(); w=t.number_per_msun.to_numpy()
        ridge=report['cmd_ridge']
        retained=T.cmd_select(pd.DataFrame(dict(phot_g_mean_mag=g,bp_rp=colour)),np.array(ridge['G']),np.array(ridge['color']),np.array(ridge['width']))
        total=float(w[hb].sum())
        metrics=dict(eta_reimers=eta,n_hb_per_represented_living_msun=total,
            model_hb_blue_box_fraction=None if total==0 else float(w[hb&in_box(g,colour)].sum()/total),
            noiseless_cmd_only_hb_retention=None if total==0 else float(w[hb&retained].sum()/total),
            model_hb_color_range=None if not hb.any() else [float(colour[hb].min()),float(colour[hb].max())],
            source_path=str(path.resolve()),source_sha256=sha256_file(path),source=imeta)
        pop.metadata['official_model']=imeta
        models.append((pop,colour,g,hb,metrics))
    out.mkdir(parents=True)
    for pop,_,_,_,m in models: pop.save(out/f"population_eta_{m['eta_reimers']:g}")
    summary=dict(n_core_stars=len(d),observed_stars_in_illustrative_blue_hb_box=observed_box,box=box,
        observed_box_stars_retained_by_original_cmd=observed_box_retained,
        cmd_cut_applied_to_observed_sample=False,source_field=str(field),source_sha256=sha256_file(field),
        distance_kpc=distance,ebv=ebv,radius_deg=radius,models=[m[-1] for m in models],
        interpretation='Illustrative sensitivity to RGB stellar winds. Eta_Reimers is NOT the globular-cluster escape rate. No eta, membership probability or age is fitted.',
        limitations=['The box is an illustrative colour/magnitude region, not a spectroscopic classification or an HB-index measurement.',
            'Core crowding/completeness, binaries, metallicity/helium/age distributions and error convolution are not calibrated here.',
            'Noiseless CMD-only retention is not a Gaia selection completeness or a mass-loss measurement.'])
    (out/'hb_comparison.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    fig,axes=plt.subplots(1,len(models),figsize=(5*len(models),6),sharex=True,sharey=True,layout='constrained')
    for ax,(pop,c,g,hb,m) in zip(np.atleast_1d(axes),models):
        ax.scatter(d.bp_rp,d.phot_g_mean_mag,s=8,color='#888888',alpha=.45,label='Core data before CMD cut')
        ax.scatter(c[~hb],g[~hb],s=9,color='#669dc8',label='Other isochrone nodes')
        ax.scatter(c[hb],g[hb],s=18,color='#d96e24',label='Helium-burning nodes (not star counts)')
        ax.add_patch(Rectangle((box['color_min'],box['g_min']),box['color_max']-box['color_min'],box['g_max']-box['g_min'],fill=False,ls='--',color='#7a5aad'))
        ax.set(xlabel='Gaia BP - RP',xlim=(-.6,2.1),ylim=(21,12),title=f"RGB wind eta = {m['eta_reimers']:g}")
        ax.grid(alpha=.15)
    np.atleast_1d(axes)[0].set_ylabel('Gaia G (mag)')
    np.atleast_1d(axes)[-1].legend(fontsize=7,loc='lower left')
    fig.suptitle('HB sensitivity at fixed age/composition/distance; illustrative box, no parameter fit',fontsize=12)
    fig.savefig(out/'hb_comparison.png',dpi=160); plt.close(fig)
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    choice=p.add_mutually_exclusive_group(required=True)
    choice.add_argument('--population',type=Path); choice.add_argument('--isochrones',type=Path,nargs='+')
    for key in ('selection','out'): p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    result=run(a.population,a.selection,a.out) if a.population else compare_isochrones(a.isochrones,a.selection,a.out)
    print('Core stars before CMD cut:',result['n_core_stars'])
