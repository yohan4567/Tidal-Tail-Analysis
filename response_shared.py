"""Share identical observation draws across quality states and weight-only populations.

No reduction of calibration states or measurement-noise draws. The node ordering,
photometry and positive-weight support must match the scalar reference exactly.
"""
import numpy as np
import pandas as pd
from response_model import science_mask


def expected_shared(particles,cases,survey,errors,T,row,P,meta,inside,edges,*,a_g,a_bp,a_rp,weights,
                    draws=16,seed=42,parallax_zero_point_mas=0.,max_pairs=10000,diagnostics=None):
    if not cases:raise ValueError('At least one population/quality case required')
    if int(draws)!=draws or draws<2 or int(max_pairs)!=max_pairs or max_pairs<1:raise ValueError('Invalid observation draws or chunk size')
    if not np.isfinite([a_g,a_bp,a_rp,parallax_zero_point_mas]).all() or min(a_g,a_bp,a_rp)<0:raise ValueError('Invalid extinction or zero point')
    required=['ra_deg','dec_deg','distance_kpc','pmra_masyr','pmdec_masyr'];weights=np.asarray(weights,float)
    if weights.shape!=(len(particles),):raise ValueError('One scalar weight per particle required')
    if not np.isfinite(particles[required]).all().all() or (particles.distance_kpc<=0).any() or not np.isfinite(weights).all() or (weights<0).any():raise ValueError('Invalid particle data or weights')
    pop=cases[0]['population'].table;basis=cases[0]['population'].metadata.get('mass_basis')
    node_table=pop.drop(columns='number_per_msun');support=pop.number_per_msun.to_numpy()>0
    for case in cases:
        other=case['population']
        if other.metadata.get('mass_basis')!=basis or not other.table.drop(columns='number_per_msun').equals(node_table):
            raise ValueError('Shared cases require identical ordered population nodes and mass basis')
        if not np.array_equal(other.table.number_per_msun.to_numpy()>0,support):raise ValueError('Positive population support differs; independent RNG streams cannot be merged')
    k=len(cases)
    if diagnostics is not None and len(diagnostics)!=k:raise ValueError('One diagnostic dictionary per case required')
    x,y=T.gnomonic(particles.ra_deg.to_numpy(),particles.dec_deg.to_numpy(),row.ra,row.dec)
    ix=np.searchsorted(edges,x,side='right')-1;iy=np.searchsorted(edges,y,side='right')-1;n=len(edges)-1
    geometry=(ix>=0)&(ix<n)&(iy>=0)&(iy<n)
    geometry &= inside[np.clip(ix,0,n-1),np.clip(iy,0,n-1)]
    a,d,a0,d0=map(np.radians,(particles.ra_deg,particles.dec_deg,row.ra,row.dec))
    sep=np.degrees(np.arccos(np.clip(np.sin(d0)*np.sin(d)+np.cos(d0)*np.cos(d)*np.cos(a-a0),-1,1)))
    geometry &= sep<=meta['radius_deg']
    for _,mx,my,mr in P.get('masks',[]):geometry &= np.hypot(x-mx,y-my)>mr
    missing=P.get('missing_tiles')
    if missing:
        pixels=T.skymath.ang2pix_nest(missing['level'],particles.ra_deg.to_numpy(),particles.dec_deg.to_numpy())
        geometry &= ~np.isin(pixels,missing['tiles'])
    if P['drift']['corrected']:
        pm,dexp,_,_=T.orbit_expectation(row,P['R'],x,y);dg=5*np.log10(dexp/row.R_sun)
    else:
        pm=(np.full(len(x),row.pmra),np.full(len(x),row.pmdec));dexp=np.full(len(x),row.R_sun);dg=None
    flow=[]
    for case in cases:
        total=case['population'].table.number_per_msun.sum()
        flow.append(dict(before_footprint=float(np.sum(weights)*total),after_footprint=float(np.sum(weights[geometry])*total),
            after_photometric_support=0.,after_gaia=0.,after_availability=0.,after_quality=0.,after_science=0.))
    realizations=np.zeros((k,draws,n,n));tracer=np.zeros((k,len(particles))) if diagnostics is not None else None
    rng=np.random.default_rng(seed);n_pop=len(pop);step=max(1,max_pairs//n_pop);geometry_ids=np.flatnonzero(geometry)
    for start in range(0,len(geometry_ids),step):
        ids=np.repeat(geometry_ids[start:start+step],n_pop);lp=np.tile(np.arange(n_pop),len(ids)//n_pop)
        dist=particles.distance_kpc.to_numpy()[ids]
        g=pop.M_G.to_numpy()[lp]+5*np.log10(dist)+10+a_g
        bp=pop.M_BP.to_numpy()[lp]+5*np.log10(dist)+10+a_bp
        rp=pop.M_RP.to_numpy()[lp]+5*np.log10(dist)+10+a_rp
        if all(c in pop for c in ('A_G','A_BP','A_RP')):
            g+=pop.A_G.to_numpy()[lp];bp+=pop.A_BP.to_numpy()[lp];rp+=pop.A_RP.to_numpy()[lp]
        reference_weight=weights[ids]*pop.number_per_msun.to_numpy()[lp]
        possible=(reference_weight>0)&(g<=min(T.G_MAX,meta['g_limit'])+8*errors.max_g_error)
        if len(P['mids'])>=3:
            lower=P['mids'].min()-.3+(0 if dg is None else dg[ids]);possible &= g>=lower-8*errors.max_g_error
        ids,lp,g,bp,rp=[v[possible] for v in (ids,lp,g,bp,rp)]
        if not len(ids):continue
        xx,yy=x[ids],y[ids]
        pair=[weights[ids]*case['population'].table.number_per_msun.to_numpy()[lp] for case in cases]
        pg=survey.query(particles.ra_deg.to_numpy()[ids],particles.dec_deg.to_numpy()[ids],g,xx,yy)
        for j in range(k):
            flow[j]['after_photometric_support']+=float(pair[j].sum());flow[j]['after_gaia']+=float(np.sum(pair[j]*pg))
        active=pg>0
        if not active.any():continue
        ids,g,bp,rp,pg,xx,yy=[v[active] for v in (ids,g,bp,rp,pg,xx,yy)];pair=[v[active] for v in pair]
        selected=[]
        for j,case in enumerate(cases):
            pa,pq=case['quality'].query(xx,yy,g,g-rp)
            flow[j]['after_availability']+=float(np.sum(pair[j]*pg*pa))
            value=pair[j]*pg*pa*pq;selected.append(value);flow[j]['after_quality']+=float(value.sum())
        support_ids=errors.donor_support(g,g-rp,xx,yy)
        for obs in range(draws):
            roots,phot_sigma=errors.sample(g,g-rp,xx,yy,rng,support=support_ids)
            delta=np.einsum('nij,nj->ni',roots,rng.normal(size=(len(ids),3)))
            mags=np.column_stack([g,bp,rp])+phot_sigma*np.clip(rng.normal(size=(len(ids),3)),-8,8)
            data=pd.DataFrame(dict(parallax=1/particles.distance_kpc.to_numpy()[ids]+parallax_zero_point_mas+delta[:,0],
                pmra=particles.pmra_masyr.to_numpy()[ids]+delta[:,1],pmdec=particles.pmdec_masyr.to_numpy()[ids]+delta[:,2],
                parallax_error=np.linalg.norm(roots[:,0,:],axis=1),pmra_error=np.linalg.norm(roots[:,1,:],axis=1),
                pmdec_error=np.linalg.norm(roots[:,2,:],axis=1),phot_g_mean_mag=mags[:,0],bp_rp=mags[:,1]-mags[:,2]))
            keep=science_mask(T,data,P,meta,(pm[0][ids],pm[1][ids]),dexp[ids],None if dg is None else dg[ids])
            for j in range(k):
                np.add.at(realizations[j,obs],(ix[ids[keep]],iy[ids[keep]]),selected[j][keep])
                if tracer is not None:np.add.at(tracer[j],ids[keep],selected[j][keep]/draws)
    outputs=[]
    for j in range(k):
        signal=realizations[j].mean(axis=0);flow[j]['after_science']=float(signal.sum())
        if diagnostics is not None:
            fit=diagnostics[j].get('fitmask',inside)[np.clip(ix,0,n-1),np.clip(iy,0,n-1)]&geometry
            diagnostics[j]['tracer_selected_counts']=tracer[j]*fit
        outputs.append((signal,realizations[j].std(axis=0,ddof=1)/np.sqrt(draws),flow[j]))
    return outputs
