"""Mass-normalized expected Gaia counts using a population and conditional selection model."""
import argparse
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import shutil
from importlib.metadata import version

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from population_model import Population
from provenance import sha256_file
from selection_bridge import load_study, DEFAULT_STUDY
from selection_model import Completeness, QualityMap, normalize_gaia, science_available, quality_without_magnitude

HERE=Path(__file__).resolve().parent
CORRELATIONS=('parallax_pmra_corr','parallax_pmdec_corr','pmra_pmdec_corr')
PHOT_ERRORS=tuple(f'phot_{band}_mean_mag_error' for band in ('g','bp','rp'))


def astrometric_covariance(df,policy='required'):
    sigma=df[['parallax_error','pmra_error','pmdec_error']].to_numpy(float)
    if not np.isfinite(sigma).all() or (sigma<0).any():
        raise ValueError('Finite nonnegative astrometric errors required')
    r=np.broadcast_to(np.eye(3),(len(df),3,3)).copy()
    if policy=='required':
        if any(k not in df for k in CORRELATIONS):
            raise ValueError('Astrometric correlation columns missing; obtain them or explicitly choose diagonal_approximation')
        for key,(i,j) in zip(CORRELATIONS,((0,1),(0,2),(1,2))):
            c=df[key].to_numpy(float)
            if not np.isfinite(c).all() or (np.abs(c)>1).any():
                raise ValueError('Invalid correlation coefficient')
            r[:,i,j]=r[:,j,i]=c
    elif policy!='diagonal_approximation':
        raise ValueError('Unknown covariance policy')
    values,vectors=np.linalg.eigh(r)
    if (values< -1e-10).any():
        raise ValueError('Astrometric correlation matrix is not positive semidefinite')
    root=vectors*np.sqrt(np.maximum(values,0))[:,None,:]
    # Root transforms standard normals into errors in mas, mas/yr, mas/yr.
    return sigma[:,:,None]*root


class ErrorLibrary:
    def __init__(self,parent,T,metadata,spec,center):
        if not metadata.get('parent_is_unfiltered',False) or not metadata.get('source_description','').strip():
            raise ValueError('Error donors require documented unfiltered-parent metadata')
        d=normalize_gaia(parent).reset_index(drop=True)
        available=science_available(d)
        d=quality_without_magnitude(T,d.loc[available]).copy()
        before_aux=len(d)
        if any(k not in d for k in PHOT_ERRORS):
            raise ValueError('Need G/BP/RP flux SNRs or explicit phot_*_mean_mag_error columns; no zero-noise substitution')
        cols=list(PHOT_ERRORS)+([] if spec.get('covariance_policy','required')=='diagonal_approximation' else list(CORRELATIONS))
        if any(k not in d for k in cols):
            raise ValueError('Error library lacks required covariance/photometry columns')
        d=d.loc[np.isfinite(d[cols].to_numpy(float)).all(axis=1)].reset_index(drop=True)
        if len(d)<int(spec.get('min_donors',20)) or (d[list(PHOT_ERRORS)]<0).any().any():
            raise ValueError('Insufficient or invalid error donors')
        self.d=d
        self.metadata=dict(n_parent=len(parent),n_science_available=int(available.sum()),n_joint_quality=before_aux,
                           n_auxiliary_error_complete=len(d),
                           caveat='Auxiliary-error missingness is not another science cut; transfer of this donor distribution must be checked.')
        self.roots=astrometric_covariance(d,spec.get('covariance_policy','required'))
        self.phot=d[list(PHOT_ERRORS)].to_numpy(float)
        self.max_g_error=float(self.phot[:,0].max())
        x,y=T.gnomonic(d.ra.to_numpy(),d.dec.to_numpy(),*center)
        self.scale=np.asarray(spec.get('feature_scales',[.5,.3,5.,5.]),float)
        if self.scale.shape!=(4,) or not np.isfinite(self.scale).all() or (self.scale<=0).any():
            raise ValueError('Four positive feature scales [G,G-RP,x,y] required')
        features=np.column_stack([d.phot_g_mean_mag,d.g_rp,x,y])/self.scale
        self.tree=cKDTree(features)
        self.radius=float(spec.get('max_feature_distance',2.))
        self.neighbours=min(int(spec.get('neighbours',20)),len(d))
        if self.radius<=0 or self.neighbours<1:
            raise ValueError('Positive donor radius/neighbour count required')

    def donor_support(self,g,c,x,y):
        features=np.column_stack([g,c,x,y])/self.scale
        distance,index=self.tree.query(features,k=self.neighbours)
        distance,index=np.asarray(distance).reshape(len(g),-1),np.asarray(index).reshape(len(g),-1)
        count=(distance<=self.radius).sum(axis=1)
        if (count==0).any():
            raise ValueError('No error donors within calibrated feature support; extend calibration or reduce analysis range')
        return index,count

    def sample(self,g,c,x,y,rng,support=None):
        index,count=self.donor_support(g,c,x,y) if support is None else support
        choice=(rng.random(len(g))*count).astype(int)
        ix=index[np.arange(len(g)),choice]
        return self.roots[ix],self.phot[ix]


def history_weights(particles,run,history=None):
    """Uniform release-time importance weights; this never evolves the progenitor potential."""
    if history is not None and run.get('dynamics_mass_history') is not None:
        raise ValueError('Cannot reweight a rate-coupled history; rerun dynamics')
    c=run['config']; duration=float(c['duration_myr'])
    age=particles.lookback_release_myr.to_numpy(float)
    if not np.isfinite(age).all() or (age<0).any() or (age>duration+1e-8).any():
        raise ValueError('Invalid release ages')
    if history is None:
        w=particles.mass_weight_msun.to_numpy(float)
        return w,float(c['reference_mdot_msun_myr']),dict(mode='constant_reference')
    needed=['lookback_myr','escaped_stellar_rate_msun_myr']
    if any(k not in history for k in needed):
        raise ValueError('History requires lookback_myr and escaped_stellar_rate_msun_myr (stellar escape only)')
    h=history.sort_values('lookback_myr')
    t,r=h[needed[0]].to_numpy(float),h[needed[1]].to_numpy(float)
    if len(t)<2 or not np.isfinite([t,r]).all() or (np.diff(t)<=0).any() or (r<0).any() or t[0]>0 or t[-1]<duration:
        raise ValueError('Finite nonnegative rate history must cover [0,duration] with increasing times')
    knots=np.unique(np.r_[0,t[(t>0)&(t<duration)],duration])
    integral=float(np.trapezoid(np.interp(knots,t,r),knots))
    if integral<=0:
        raise ValueError('Rate history has zero escaped mass')
    # Official sampler in this project releases uniformly; T/N accounts for both arms.
    w=duration/len(particles)*np.interp(age,t,r)
    return w,integral/duration,dict(mode='time_reweighted_fixed_dynamics',integrated_escaped_mass_msun=integral,
                                   sampled_mass_msun=float(w.sum()),
                                   caveat='No feedback of changed escape history on M(t), radius or particle release kinematics.')


def science_mask(T,df,P,meta,pm_expected,d_expected,dg):
    """Only science/download cuts: availability and quality have already been weighted."""
    keep=(df.phot_g_mean_mag<T.G_MAX)&(df.phot_g_mean_mag<meta['g_limit'])
    keep &= (df.pmra>=meta['pm0'][0]-meta['pm_box'])&(df.pmra<=meta['pm0'][0]+meta['pm_box'])
    keep &= (df.pmdec>=meta['pm0'][1]-meta['pm_box'])&(df.pmdec<=meta['pm0'][1]+meta['pm_box'])
    keep &= T.plx_select(df,d_expected)&T.pm_select(df,pm_expected,d_expected,P['sig0'])
    keep &= T.cmd_select(df,P['mids'],P['med'],P['wid'],dg)
    return np.asarray(keep,bool)


def expected_template(particles,population,quality,survey,errors,T,row,P,meta,inside,edges,*,
                      a_g,a_bp,a_rp,weights,draws=16,seed=42,parallax_zero_point_mas=0.,max_pairs=10000,diagnostics=None):
    if int(draws)!=draws or draws<2 or int(max_pairs)!=max_pairs or max_pairs<1:
        raise ValueError('At least 2 observation draws and positive integer chunk size required')
    if not np.isfinite([a_g,a_bp,a_rp,parallax_zero_point_mas]).all() or min(a_g,a_bp,a_rp)<0:
        raise ValueError('Finite extinction/zero point and nonnegative extinction required')
    required=['ra_deg','dec_deg','distance_kpc','pmra_masyr','pmdec_masyr']
    if not np.isfinite(particles[required]).all().all() or (particles.distance_kpc<=0).any() or not np.isfinite(weights).all() or (weights<0).any():
        raise ValueError('Invalid dynamical particle data/weights')
    x,y=T.gnomonic(particles.ra_deg.to_numpy(),particles.dec_deg.to_numpy(),row.ra,row.dec)
    ix=np.searchsorted(edges,x,side='right')-1; iy=np.searchsorted(edges,y,side='right')-1
    n=len(edges)-1
    geometry=(ix>=0)&(ix<n)&(iy>=0)&(iy<n)
    geometry &= inside[np.clip(ix,0,n-1),np.clip(iy,0,n-1)]
    a,d,a0,d0=map(np.radians,(particles.ra_deg,particles.dec_deg,row.ra,row.dec))
    sep=np.degrees(np.arccos(np.clip(np.sin(d0)*np.sin(d)+np.cos(d0)*np.cos(d)*np.cos(a-a0),-1,1)))
    geometry &= sep<=meta['radius_deg']
    for _,mx,my,mr in P.get('masks',[]):
        geometry &= np.hypot(x-mx,y-my)>mr
    missing=P.get('missing_tiles')
    if missing:
        pixels=T.skymath.ang2pix_nest(missing['level'],particles.ra_deg.to_numpy(),particles.dec_deg.to_numpy())
        geometry &= ~np.isin(pixels,missing['tiles'])
    if P['drift']['corrected']:
        pm,dexp,_,_=T.orbit_expectation(row,P['R'],x,y)
        dg=5*np.log10(dexp/row.R_sun)
    else:
        pm=(np.full(len(x),row.pmra),np.full(len(x),row.pmdec)); dexp=np.full(len(x),row.R_sun); dg=None
    pop=population.table
    total_stars=float(np.sum(weights)*pop.number_per_msun.sum())
    cutflow=dict(before_footprint=total_stars,after_footprint=float(np.sum(weights[geometry])*pop.number_per_msun.sum()),
                 after_photometric_support=0.,after_gaia=0.,after_availability=0.,after_quality=0.,after_science=0.)
    rng=np.random.default_rng(seed)
    realizations=np.zeros((draws,n,n))
    tracer_counts=np.zeros(len(particles)) if diagnostics is not None else None
    n_pop=len(pop)
    step=max(1,max_pairs//n_pop)
    geometry_ids=np.flatnonzero(geometry)
    for start in range(0,len(geometry_ids),step):
        ids=np.repeat(geometry_ids[start:start+step],n_pop)
        lp=np.tile(np.arange(n_pop),len(ids)//n_pop)
        dist=particles.distance_kpc.to_numpy()[ids]
        g=pop.M_G.to_numpy()[lp]+5*np.log10(dist)+10+a_g
        bp=pop.M_BP.to_numpy()[lp]+5*np.log10(dist)+10+a_bp
        rp=pop.M_RP.to_numpy()[lp]+5*np.log10(dist)+10+a_rp
        if all(k in pop for k in ('A_G','A_BP','A_RP')):
            g+=pop.A_G.to_numpy()[lp]; bp+=pop.A_BP.to_numpy()[lp]; rp+=pop.A_RP.to_numpy()[lp]
        pair_weight=weights[ids]*pop.number_per_msun.to_numpy()[lp]
        # Explicit bounded-Gaussian kernel: photometric deviates are clipped at 8 sigma.
        # Nodes unable to reach the upper scientific G cut have exactly zero response
        # under this declared kernel; their mass remains in the population denominator.
        possible=(pair_weight>0)&(g<=min(T.G_MAX,meta['g_limit'])+8*errors.max_g_error)
        if len(P['mids'])>=3:
            lower=P['mids'].min()-.3+(0 if dg is None else dg[ids])
            possible &= g>=lower-8*errors.max_g_error
        ids,lp,g,bp,rp,pair_weight=[v[possible] for v in (ids,lp,g,bp,rp,pair_weight)]
        if not len(ids):
            continue
        xx,yy=x[ids],y[ids]
        cutflow['after_photometric_support']+=float(pair_weight.sum())
        pg=survey.query(particles.ra_deg.to_numpy()[ids],particles.dec_deg.to_numpy()[ids],g,xx,yy)
        cutflow['after_gaia']+=float(np.sum(pair_weight*pg))
        active=pg>0
        if not active.any():
            continue
        ids,g,bp,rp,pair_weight,pg,xx,yy=[v[active] for v in (ids,g,bp,rp,pair_weight,pg,xx,yy)]
        pa,pq=quality.query(xx,yy,g,g-rp)
        cutflow['after_availability']+=float(np.sum(pair_weight*pg*pa))
        selected_weight=pair_weight*pg*pa*pq
        cutflow['after_quality']+=float(selected_weight.sum())
        donor_support=errors.donor_support(g,g-rp,xx,yy)
        for obs in range(draws):
            roots,phot_sigma=errors.sample(g,g-rp,xx,yy,rng,support=donor_support)
            delta=np.einsum('nij,nj->ni',roots,rng.normal(size=(len(ids),3)))
            mags=np.column_stack([g,bp,rp])+phot_sigma*np.clip(rng.normal(size=(len(ids),3)),-8,8)
            data=pd.DataFrame(dict(parallax=1/particles.distance_kpc.to_numpy()[ids]+parallax_zero_point_mas+delta[:,0],
                                   pmra=particles.pmra_masyr.to_numpy()[ids]+delta[:,1],
                                   pmdec=particles.pmdec_masyr.to_numpy()[ids]+delta[:,2],
                                   parallax_error=np.linalg.norm(roots[:,0,:],axis=1),
                                   pmra_error=np.linalg.norm(roots[:,1,:],axis=1),pmdec_error=np.linalg.norm(roots[:,2,:],axis=1),
                                   phot_g_mean_mag=mags[:,0],bp_rp=mags[:,1]-mags[:,2]))
            keep=science_mask(T,data,P,meta,(pm[0][ids],pm[1][ids]),dexp[ids],None if dg is None else dg[ids])
            np.add.at(realizations[obs],(ix[ids[keep]],iy[ids[keep]]),selected_weight[keep])
            if tracer_counts is not None:
                np.add.at(tracer_counts,ids[keep],selected_weight[keep]/draws)
    signal=realizations.mean(axis=0)
    cutflow['after_science']=float(signal.sum())
    if diagnostics is not None:
        in_fit=diagnostics.get('fitmask',inside)[np.clip(ix,0,n-1),np.clip(iy,0,n-1)] & geometry
        diagnostics['tracer_selected_counts']=tracer_counts*in_fit
    return signal,realizations.std(axis=0,ddof=1)/np.sqrt(draws),cutflow


def run(selection,simulation,config,out,study=DEFAULT_STUDY,*,backend="shared",_provider=None):
    if backend not in ("shared","reference"):raise ValueError("Unknown observation backend")
    selection,simulation,config,out=map(Path,(selection,simulation,config,out))
    spec=json.loads(config.read_text(encoding='utf-8')); base=config.resolve().parent
    if 'parallax_zero_point_mas' not in spec or not spec.get('zero_point_description','').strip():
        raise ValueError('Explicit parallax_zero_point_mas and zero_point_description required, including a declared zero-correction baseline')
    report=json.loads((selection/'selection.json').read_text(encoding='utf-8'))
    runmeta=json.loads((simulation/'run.json').read_text(encoding='utf-8'))
    if report['cluster']!=runmeta['config']['name']:
        raise ValueError('Simulation/observation cluster mismatch')
    T=load_study(study,report['legacy'])
    selection_hash=sha256_file(Path(T.HERE)/'step3_tail_search.py')
    if selection_hash!=report['selection_file_sha256']:
        raise ValueError('Selection code changed; regenerate observed selection')
    if 'cluster_selection_parameters' not in report:
        raise ValueError('Regenerate selection_bridge output with frozen selection parameters')
    row=pd.Series(report['cluster_selection_parameters'],name=report['cluster'])
    ridge=report['cmd_ridge']
    P=dict(R=report['field_radius_deg'],sig0=report['sig0_kms'],drift=report['drift'],
           mids=np.array(ridge['G']),med=np.array(ridge['color']),wid=np.array(ridge['width']),
           masks=report['neighbour_masks'],missing_tiles=report.get('missing_tiles'))
    pop=Population.load(base/spec['population'])
    if spec.get('mass_basis')!=pop.metadata['mass_basis']:
        raise ValueError('Response and population mass_basis must match explicitly')
    quality_path=base/spec['quality']
    quality=QualityMap.load(quality_path)
    if quality.metadata['selection_sha256']!=selection_hash or not np.allclose([row.ra,row.dec],[quality.metadata['center_ra_deg'],quality.metadata['center_dec_deg']],atol=1e-8,rtol=0):
        raise ValueError('Quality calibration selection/sky-centre mismatch')
    survey=Completeness(spec['completeness'],base)
    error_spec=spec['errors']; error_path=base/error_spec['path']
    error_meta_path=base/error_spec['metadata']
    error_meta=json.loads(error_meta_path.read_text(encoding='utf-8'))
    errors=ErrorLibrary(pd.read_csv(error_path,low_memory=False),T,error_meta,error_spec,(row.ra,row.dec))
    particles=pd.read_csv(simulation/'particles.csv')
    history_path=base/spec['rate_history'] if spec.get('rate_history') else None
    dynamics_history=runmeta.get('dynamics_mass_history')
    if dynamics_history is not None:
        if history_path is not None:
            raise ValueError('An evolving-mass simulation cannot be reweighted to a different escape history; rerun its dynamics')
        if spec['mass_basis']!=dynamics_history['population_mass_basis']:
            raise ValueError('Population mass basis differs from the prescribed dynamical mass budget')
    if history_path is not None and (spec.get('rate_history_mass_basis')!=spec['mass_basis'] or not spec.get('rate_history_source','').strip()):
        raise ValueError('Rate history needs a matching rate_history_mass_basis and documented rate_history_source')
    if history_path is not None and runmeta.get('release_time_sampling') not in (None,'uniform_lookback_both_arms'):
        raise ValueError('Importance reweighting requires uniform lookback-time sampling')
    if history_path is not None and runmeta.get('release_time_sampling') is None and not spec.get('assume_uniform_release_for_legacy',False):
        raise ValueError('Legacy simulation needs explicit assume_uniform_release_for_legacy for time reweighting')
    weights,reference_rate,history_meta=history_weights(particles,runmeta,None if history_path is None else pd.read_csv(history_path))
    from joint_grid import window_weights
    window=spec.get('lookback_window_myr',runmeta['config']['duration_myr'])
    weights,window_meta=window_weights(particles,runmeta,window,weights)
    with np.load(selection/'maps.npz',allow_pickle=False) as d:
        maps={k:d[k] for k in ('counts','background','inside','edges','fitmask')}
    identity=hashlib.sha256()
    for name,array in sorted(maps.items()):
        identity.update(name.encode()); identity.update(str(array.shape).encode()); identity.update(str(array.dtype).encode()); identity.update(array.tobytes())
    identity.update(selection_hash.encode()); identity.update(report['source_field_sha256'].encode())
    identity.update(json.dumps({k:report[k] for k in ('source_metadata','cluster_selection_parameters','cmd_ridge','drift')},sort_keys=True).encode())
    observation_identity=identity.hexdigest()
    inputs=[config,selection/'selection.json',selection/'maps.npz',simulation/'particles.csv',simulation/'run.json',
            base/spec['population']/'population.csv',base/spec['population']/'population.json',
            quality_path/'quality.npz',quality_path/'quality.json',error_path,error_meta_path]
    if history_path is not None:
        inputs.append(history_path)
    input_hashes={str(f.resolve()):sha256_file(f) for f in inputs}
    source_hashes={name:sha256_file(HERE/name) for name in ('response_model.py','selection_model.py','population_model.py','joint_grid.py')}
    versions={name:version(name) for name in ('numpy','scipy','pandas','astropy','astropy-healpix','h5py')}
    if backend=='shared':source_hashes['response_shared.py']=sha256_file(HERE/'response_shared.py')
    if _provider is not None:source_hashes['response_group.py']=sha256_file(HERE/'response_group.py')
    cache_key=hashlib.sha256(json.dumps(dict(inputs=input_hashes,source=source_hashes,versions=versions,completeness=survey.metadata,backend=backend),sort_keys=True).encode()).hexdigest()
    cache=(base/spec['cache_dir']/cache_key).resolve() if spec.get('cache_dir') else None
    if cache is not None and (cache==out.resolve() or out.resolve() in cache.parents or cache in out.resolve().parents):
        raise ValueError('Response cache and output directories must be disjoint')
    if cache is not None and (cache/'response.json').exists():
        cached=json.loads((cache/'response.json').read_text(encoding='utf-8'))
        if cached['template_sha256']!=sha256_file(cache/'calibrated_counts.npz') or cached['maps_sha256']!=sha256_file(cache/'response_maps.npz'):
            raise ValueError('Response cache integrity check failed')
        shutil.copytree(cache,out)
        cached.update(cache_hit=True,cache_key=cache_key)
        (out/'response.json').write_text(json.dumps(cached,indent=2),encoding='utf-8')
        print(f'Reused response cache {cache_key[:12]}',flush=True)
        return cached
    count=int(spec.get('calibration_draws',1)); draws=int(spec.get('observation_draws',16)); seed=int(spec.get('seed',42))
    if count<1:
        raise ValueError('At least one calibration draw required')
    cases=[dict(population=pop,quality=quality if count==1 else QualityMap.load(quality_path,seed+100000+k)) for k in range(count)]
    diagnostics=[dict(fitmask=maps['fitmask']) for _ in cases]
    options=dict(weights=weights,draws=draws,seed=seed,
        a_g=spec['extinction']['a_g'],a_bp=spec['extinction']['a_bp'],a_rp=spec['extinction']['a_rp'],
        parallax_zero_point_mas=float(spec['parallax_zero_point_mas']),max_pairs=int(spec.get('max_pairs',10000)))
    args=(particles,cases,survey,errors,T,row,P,report['source_metadata'],maps['inside'],maps['edges'])
    if backend=='shared':
        from response_shared import expected_shared
        results=(_provider or expected_shared)(*args,diagnostics=diagnostics,**options)
    else:
        results=[expected_template(particles,c['population'],c['quality'],survey,errors,T,row,P,report['source_metadata'],
            maps['inside'],maps['edges'],diagnostics=d,**options) for c,d in zip(cases,diagnostics)]
    signals=[]; mc_errors=[]; cutflows=[]; tracer_counts=[]
    for k,((signal,mc_error,flow),diagnostic) in enumerate(zip(results,diagnostics)):
        tracer_counts.append(diagnostic['tracer_selected_counts'])
        signals.append(signal); mc_errors.append(mc_error); cutflows.append(flow)
        print(f'Observation response {k+1}/{count}: {signal.sum():.3f} expected selected stars',flush=True)
    signals=np.asarray(signals)
    synthetic=any(bool(m.get('synthetic_calibration',False)) for m in (spec,pop.metadata,quality.metadata,survey.metadata,error_meta))
    provenance=dict(schema_version=2,cluster=report['cluster'],mass_basis=spec['mass_basis'],reference_mdot=reference_rate,
        response_sampling_identity=hashlib.sha256(json.dumps({k:v for k,v in spec.items() if k not in ('lookback_window_myr','cache_dir')},sort_keys=True).encode()).hexdigest(),
        dynamics_sampling_identity=str(runmeta['config']['seed']),
        dynamics_family_identity=hashlib.sha256(json.dumps({k:v for k,v in runmeta['config'].items() if k not in ('seed','duration_myr','reference_mdot_msun_myr','n_particles')},sort_keys=True).encode()).hexdigest(),
        model_version=runmeta.get('model_version','legacy_unversioned'),lookback_window=window_meta,
        dynamics_mass_history=dynamics_history,amplitude_policy='fixed_reference_rate' if dynamics_history is not None else 'free_amplitude_fixed_dynamics',
        scientific_status=spec.get('scientific_status','conditional_model'),scenario=spec.get('scenario'),
        rate_scope=pop.metadata.get('rate_scope','Mass represented by supplied population only'),
        extinction=spec['extinction'],extinction_description=spec.get('extinction_description','Declared constant extinction'),
        synthetic_calibration=synthetic,population=pop.metadata,quality=quality.metadata,completeness=survey.metadata,
        error_source=error_meta,error_library_summary=errors.metadata,covariance_policy=error_spec.get('covariance_policy','required'),
        parallax_zero_point_mas=spec['parallax_zero_point_mas'],zero_point_description=spec['zero_point_description'],
        cutflow_expected_counts=cutflows,rate_history=history_meta,
        observation_draws=draws,calibration_draws=count,seed=seed,
        input_sha256=input_hashes,source_sha256=source_hashes,versions=versions,observation_identity=observation_identity,
        cache_key=cache_key,cache_hit=False,observation_backend=("shared_group_v1" if _provider is not None else "shared_quality_v1") if backend=="shared" else "reference_v1",
        limitations=['Quality/availability calibration transfer and its covariate sufficiency require validation.',
                     'Completeness/quality lookup uses noiseless apparent photometry as a small-noise approximation.',
                     'Photometric Gaussian draws are clipped at 8 sigma, with no cross-band error covariance.',
                     'No spatially correlated astrometric systematics; zero point supplied explicitly.',
                     'Only availability/quality Beta-cell uncertainty is sampled automatically; other calibration variants require separate configurations.',
                     'Background shape remains the original fitted map; background mismatch must be checked separately.',
                     'Coupled histories use fixed rate amplitude; fixed-dynamics reweighting cannot change progenitor mass/size.',
                     'Finite tracer and observation sampling must pass convergence checks before physical inference.'])
    out.mkdir(parents=True,exist_ok=False)
    average=signals.mean(axis=0); fitmask=maps['fitmask']
    selection_description=json.dumps(dict(response_model='hybrid_v2',mass_basis=spec['mass_basis'],synthetic_calibration=synthetic,
                                         response_sampling_identity=provenance['response_sampling_identity'],dynamics_sampling_identity=provenance['dynamics_sampling_identity'],
                                         dynamics_family_identity=provenance['dynamics_family_identity'],
                                         model_version=provenance['model_version'],lookback_window=window_meta,
                                         dynamics_mass_history=dynamics_history,amplitude_policy=provenance['amplitude_policy'],
                                         scientific_status=provenance['scientific_status'],rate_scope=provenance['rate_scope'],scenario=provenance['scenario'],
                                         rate_history=history_meta,
                                         observation_identity=observation_identity,
                                         selection_sha256=selection_hash,source_description=spec.get('source_description',''),limitations=provenance['limitations']))
    np.savez_compressed(out/'calibrated_counts.npz',counts=maps['counts'][fitmask],background=maps['background'][fitmask],
                        signal=average[fitmask],signal_draws=signals[:,fitmask],reference_mdot=reference_rate,
                        synthetic_calibration=synthetic,selection_description=selection_description)
    np.savez_compressed(out/'response_maps.npz',signal=average,signal_draws=signals,mc_standard_error=np.mean(mc_errors,axis=0),**maps)
    ages=particles.lookback_release_myr.to_numpy(float)
    age_edges=np.linspace(0,runmeta['config']['duration_myr'],21)
    selected=np.mean(tracer_counts,axis=0)
    age_counts=np.histogram(ages,bins=age_edges,weights=selected)[0]
    pd.DataFrame(dict(lookback_start_myr=age_edges[:-1],lookback_end_myr=age_edges[1:],
        expected_selected_stars=age_counts,cumulative_older_stars=np.cumsum(age_counts[::-1])[::-1])).to_csv(out/'age_contributions.csv',index=False)
    tx,ty=T.gnomonic(particles.ra_deg.to_numpy(),particles.dec_deg.to_numpy(),row.ra,row.dec)
    nx=len(maps['edges'])-1
    pixx=np.searchsorted(maps['edges'],tx,side='right')-1; pixy=np.searchsorted(maps['edges'],ty,side='right')-1
    lookup=np.full(maps['fitmask'].shape,-1,int); lookup[maps['fitmask']]=np.arange(maps['fitmask'].sum())
    valid=(pixx>=0)&(pixx<nx)&(pixy>=0)&(pixy<nx)
    fit_pixel=np.full(len(particles),-1,int); fit_pixel[valid]=lookup[pixx[valid],pixy[valid]]
    np.savez_compressed(out/'tracer_response.npz',lookback_myr=ages,expected_selected_stars=selected,
        mass_weight_msun=weights,calibration_counts=np.asarray(tracer_counts),fit_pixel_index=fit_pixel)
    provenance['expected_stars_in_fit_region']=float(average[fitmask].sum())
    provenance['template_sha256']=sha256_file(out/'calibrated_counts.npz')
    provenance['maps_sha256']=sha256_file(out/'response_maps.npz')
    (out/'response.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    os.environ.setdefault('MPLCONFIGDIR',str(HERE/'outputs/.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(14,4.5),constrained_layout=True)
    extent=[maps['edges'][0],maps['edges'][-1],maps['edges'][0],maps['edges'][-1]]
    for ax,image,title in zip(axes,[maps['counts'],average,np.mean(mc_errors,axis=0)],
                             ['Observed candidate counts','Reference-rate expected signal','Observation Monte Carlo standard error']):
        im=ax.imshow(np.where(maps['inside'],image,np.nan).T,origin='lower',extent=extent)
        ax.set(title=title,xlabel='Gnomonic east (deg)',ylabel='Gnomonic north (deg)')
        fig.colorbar(im,ax=ax,label='Stars / pixel')
    fig.suptitle(f"{report['cluster']} | {'SYNTHETIC calibration test' if synthetic else 'conditional calibrated model'}")
    fig.savefig(out/'response.png',dpi=150); plt.close(fig)
    if cache is not None and not cache.exists():
        cache.parent.mkdir(parents=True,exist_ok=True)
        shutil.copytree(out,cache)
    return provenance


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('selection','run','config','out'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--study',type=Path,default=DEFAULT_STUDY)
    p.add_argument('--backend',choices=['shared','reference'],default='shared')
    a=p.parse_args()
    result=run(a.selection,a.run,a.config,a.out,a.study,backend=a.backend)
    print(json.dumps(dict(expected_stars_in_fit_region=result['expected_stars_in_fit_region'],mass_basis=result['mass_basis'],synthetic=result['synthetic_calibration']),indent=2))


if __name__=='__main__':
    main()
