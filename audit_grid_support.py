"""Audit every existing simulation using the frozen response's support calculations.

This diagnostic omits ONLY the observation-draw loop. Its zero signal arrays are
discarded and are never exported as a response or used in a likelihood.
"""
import argparse,ast,inspect,json,os,sys
from pathlib import Path


def support_function(function):
    tree=ast.parse(inspect.getsource(function));removed=[]
    class OmitDraws(ast.NodeTransformer):
        def visit_For(self,node):
            if isinstance(node.target,ast.Name) and node.target.id=='obs':
                if ast.unparse(node.iter)!='range(draws)':raise ValueError('Unrecognized observation loop')
                removed.append(node);return ast.copy_location(ast.Pass(),node)
            return self.generic_visit(node)
    tree=OmitDraws().visit(tree)
    if len(removed)!=1:raise ValueError('Expected exactly one frozen observation-draw loop')
    ast.fix_missing_locations(tree);namespace={}
    exec(compile(tree,'<support-only audit; no response output>','exec'),function.__globals__,namespace)
    return namespace[function.__name__]


def run(job,grid,out):
    from provenance import sha256_file
    j=json.loads(Path(job).read_text(encoding='utf-8-sig'))
    os.environ['GC_TIDAL_STUDY']=j['study'];sys.path.insert(0,j['code'])
    import numpy as np,pandas as pd
    from response_model import ErrorLibrary,Completeness,Population,QualityMap,load_study
    from response_shared import expected_shared
    from empirical_population import prepare
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    prior_path=out/'population/population_prior.json'
    prior=json.loads(prior_path.read_text()) if prior_path.exists() else prepare(j['cluster'],j['isochrone'],j['reference'],j['html'],prior_path.parent,offset=0.,order=5)
    config=Path(prior['response_configs'][2]);spec=json.loads(config.read_text());base=config.parent
    selection=Path(j['selection']);report=json.loads((selection/'selection.json').read_text())
    T=load_study(Path(j['study']),report['legacy'])
    if sha256_file(Path(T.HERE)/'step3_tail_search.py')!=report['selection_file_sha256']:raise ValueError('Selection changed')
    row=pd.Series(report['cluster_selection_parameters'],name=report['cluster']);ridge=report['cmd_ridge']
    P=dict(R=report['field_radius_deg'],sig0=report['sig0_kms'],drift=report['drift'],mids=np.array(ridge['G']),
           med=np.array(ridge['color']),wid=np.array(ridge['width']),masks=report['neighbour_masks'],missing_tiles=report.get('missing_tiles'))
    with np.load(selection/'maps.npz') as f:inside,edges=f['inside'],f['edges']
    quality=QualityMap.load(base/spec['quality']);pop=Population.load(base/spec['population'])
    survey=Completeness(spec['completeness'],base)
    error_meta=json.loads((base/spec['errors']['metadata']).read_text())
    errors=ErrorLibrary(pd.read_csv(base/spec['errors']['path'],low_memory=False),T,error_meta,spec['errors'],(row.ra,row.dec))
    evaluate=support_function(expected_shared);plan=json.loads((Path(grid)/'plan.json').read_text());results=[]
    cases=[dict(population=pop,quality=quality)]
    for rate in plan['rates']:
        for seed in plan['seeds']:
            sim=Path(grid)/f'rate_{rate:g}_seed_{seed}'/'simulation'
            particles=pd.read_csv(sim/'particles.csv')
            if (particles.mass_weight_msun<=0).any():raise ValueError('Positive full-window tracer weights required for subset coverage')
            try:
                evaluate(particles,cases,survey,errors,T,row,P,report['source_metadata'],inside,edges,
                    a_g=spec['extinction']['a_g'],a_bp=spec['extinction']['a_bp'],a_rp=spec['extinction']['a_rp'],
                    weights=particles.mass_weight_msun.to_numpy(),draws=2)
                record=dict(rate=rate,seed=seed,passed=True,particles_sha256=sha256_file(sim/'particles.csv'))
            except ValueError as exc:record=dict(rate=rate,seed=seed,passed=False,error=str(exc))
            results.append(record)
            report_out=dict(stage='all_existing_grid_support_passed' if all(r['passed'] for r in results) else 'support_failed',
                checked=len(results),expected=len(plan['rates'])*len(plan['seeds']),results=results,
                complete=len(results)==len(plan['rates'])*len(plan['seeds']),
                method='Frozen expected_shared with observation-draw loop omitted; all geometry, photometry, completeness, quality and donor-support operations retained; arrays discarded.',
                response_shared_sha256=sha256_file(Path(j['code'])/'response_shared.py'),job_sha256=sha256_file(job),
                physical_rate_measurement_validated=False)
            (out/'audit.json').write_text(json.dumps(report_out,indent=2))
            print(rate,seed,record['passed'],flush=True)
    if not all(r['passed'] for r in results):raise ValueError('Some full-grid simulations still lack calibration support')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('job','grid','out'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    run(a.job,a.grid,a.out)
