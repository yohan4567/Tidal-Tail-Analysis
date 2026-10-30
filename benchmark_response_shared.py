"""Benchmark shared observations against the reference on existing real inputs.

Writes diagnostics only; never changes a production run or its frozen plan.
"""
import argparse
import json
from pathlib import Path
from time import perf_counter
import numpy as np
import pandas as pd
from response_model import (Population, QualityMap, Completeness, ErrorLibrary,
                            expected_template, load_study, DEFAULT_STUDY, history_weights)
from response_shared import expected_shared
from joint_grid import window_weights
from provenance import sha256_file


def run(selection, simulation, configs, output, particle_count=128, study=DEFAULT_STUDY):
    if particle_count<1:raise ValueError('Positive benchmark particle count required')
    selection,simulation,output=map(Path,(selection,simulation,output))
    specs=[json.loads(Path(c).read_text(encoding='utf-8')) for c in configs]
    bases=[Path(c).resolve().parent for c in configs]
    common=('quality','mass_basis','completeness','errors','extinction','observation_draws',
            'calibration_draws','seed','parallax_zero_point_mas','lookback_window_myr','rate_history')
    for spec in specs[1:]:
        if any(spec.get(k)!=specs[0].get(k) for k in common):
            raise ValueError('Configs differ in observation settings; cannot share')
    spec=specs[0];base=bases[0]
    report=json.loads((selection/'selection.json').read_text(encoding='utf-8'))
    runmeta=json.loads((simulation/'run.json').read_text(encoding='utf-8'))
    T=load_study(study,report['legacy'])
    if sha256_file(Path(T.HERE)/'step3_tail_search.py')!=report['selection_file_sha256']:
        raise ValueError('Selection code changed')
    row=pd.Series(report['cluster_selection_parameters'],name=report['cluster']);ridge=report['cmd_ridge']
    P=dict(R=report['field_radius_deg'],sig0=report['sig0_kms'],drift=report['drift'],
           mids=np.array(ridge['G']),med=np.array(ridge['color']),wid=np.array(ridge['width']),
           masks=report['neighbour_masks'],missing_tiles=report.get('missing_tiles'))
    errors=ErrorLibrary(pd.read_csv(base/spec['errors']['path'],low_memory=False),T,
        json.loads((base/spec['errors']['metadata']).read_text(encoding='utf-8')),spec['errors'],(row.ra,row.dec))
    survey=Completeness(spec['completeness'],base)
    particles=pd.read_csv(simulation/'particles.csv')
    weights,_,_=history_weights(particles,runmeta)
    weights,_=window_weights(particles,runmeta,spec.get('lookback_window_myr',runmeta['config']['duration_myr']),weights)
    if spec.get('rate_history'):raise ValueError('Benchmark currently requires simulation-native mass history')
    maps=np.load(selection/'maps.npz');edges=maps['edges'];inside=maps['inside']
    x,y=T.gnomonic(particles.ra_deg.to_numpy(),particles.dec_deg.to_numpy(),row.ra,row.dec)
    ix=np.searchsorted(edges,x,side='right')-1;iy=np.searchsorted(edges,y,side='right')-1;n=len(edges)-1
    eligible=np.flatnonzero((ix>=0)&(ix<n)&(iy>=0)&(iy<n)&(weights>0))
    if not len(eligible):raise ValueError('No particles in map')
    ids=eligible[np.linspace(0,len(eligible)-1,min(particle_count,len(eligible)),dtype=int)]
    particles=particles.iloc[ids].reset_index(drop=True);weights=weights[ids]
    cases=[]
    for s,b in zip(specs,bases):
        pop=Population.load(b/s['population'])
        for k in range(s['calibration_draws']):
            quality=QualityMap.load(b/s['quality'],None if s['calibration_draws']==1 else s['seed']+100000+k)
            cases.append(dict(population=pop,quality=quality))
    opts=dict(**spec['extinction'],weights=weights,draws=spec['observation_draws'],seed=spec['seed'],
              parallax_zero_point_mas=spec['parallax_zero_point_mas'])
    diags=[dict(fitmask=maps['fitmask']) for _ in cases];start=perf_counter()
    reference=[expected_template(particles,c['population'],c['quality'],survey,errors,T,row,P,report['source_metadata'],inside,edges,diagnostics=d,**opts) for c,d in zip(cases,diags)]
    reference_seconds=perf_counter()-start
    shared_diags=[dict(fitmask=maps['fitmask']) for _ in cases];start=perf_counter()
    shared=expected_shared(particles,cases,survey,errors,T,row,P,report['source_metadata'],inside,edges,diagnostics=shared_diags,**opts)
    shared_seconds=perf_counter()-start
    for a,b,da,db in zip(reference,shared,diags,shared_diags):
        np.testing.assert_array_equal(a[0],b[0]);np.testing.assert_array_equal(a[1],b[1])
        if a[2]!=b[2]:raise AssertionError('Cutflow differs')
        np.testing.assert_array_equal(da['tracer_selected_counts'],db['tracer_selected_counts'])
    result=dict(particles=len(particles),cases=len(cases),draws=spec['observation_draws'],bitwise_equal=True,
        reference_seconds=reference_seconds,shared_seconds=shared_seconds,speedup=reference_seconds/shared_seconds,
        signal_counts=[float(a[0].sum()) for a in reference],
        source_sha256={p:sha256_file(p) for p in ('response_model.py','response_shared.py','benchmark_response_shared.py')},
        limitation='Subset benchmark, not a full-grid scientific validation or end-to-end speedup')
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2));return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--selection',required=True);p.add_argument('--simulation',required=True)
    p.add_argument('--configs',nargs='+',required=True);p.add_argument('--output',required=True);p.add_argument('--particles',type=int,default=128)
    a=p.parse_args();run(a.selection,a.simulation,a.configs,a.output,a.particles)
