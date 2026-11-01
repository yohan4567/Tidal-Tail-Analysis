"""Explicitly separate measured MF constraints from ESO's indirect model input."""
import json,re
from pathlib import Path
from empirical_population import prepare as measured_prepare,measured_population,normal_nodes,_Table
from provenance import sha256_file

def indirect_constraint(html,cluster):
    if cluster!='ESO452-SC11':raise ValueError('No reviewed indirect policy for '+cluster)
    parser=_Table();parser.feed(Path(html).read_text(encoding='utf-8'))
    rows=[r for r in parser.rows if len(r)>20 and re.search(r'\bESO\s*452[\s-]*SC11\b',r[0][0])]
    if len(rows)!=1 or not rows[0][20][1]:raise ValueError('Expected explicitly indirect catalogue MF row')
    values=re.findall(r'[+-]?\d+(?:\.\d+)?',rows[0][20][0].replace('\u2212','-'))
    if len(values)!=2:raise ValueError('Indirect slope/error unavailable')
    mean,sigma=map(float,values)
    if sigma<=0:raise ValueError('Positive estimator uncertainty required')
    return dict(cluster=cluster,slope_mean=mean,slope_sigma=sigma,indirect_estimate=True,
        measured_mass_range_msun=None,constraint_kind='indirect_estimator_model_assumption',
        convention='dN/dm proportional to m**slope',source_sha256=sha256_file(html),
        source='https://people.smp.uq.edu.au/HolgerBaumgardt/globular/parameter.html',paper='https://arxiv.org/abs/2303.01636',
        scope='Not a direct progenitor or escape MF measurement. Gaussian interpretation of the estimator is a conditional modelling assumption; absolute-rate gate remains closed.')

def prepare(cluster,isochrone,reference_config,html,out,*,offset=0.,order=5,allow_indirect_model=False):
    if cluster!='ESO452-SC11':
        if allow_indirect_model:raise ValueError('Indirect override is scoped to ESO452-SC11 only')
        return measured_prepare(cluster,isochrone,reference_config,html,out,offset=offset,order=order)
    if not allow_indirect_model:raise ValueError('ESO requires explicit indirect-model opt-in')
    out=Path(out);out.mkdir(parents=True,exist_ok=False);constraint=indirect_constraint(html,cluster)
    imeta=json.loads(Path(isochrone).with_suffix('.json').read_text())
    if imeta['sha256']!=sha256_file(isochrone):raise ValueError('Isochrone checksum mismatch')
    cfg=Path(reference_config).resolve();spec=json.loads(cfg.read_text());root=cfg.parent
    for key in ('quality','cache_dir','rate_history'):
        if spec.get(key):spec[key]=str((root/spec[key]).resolve())
    for key in ('path','metadata'):spec['errors'][key]=str((root/spec['errors'][key]).resolve())
    for key in ('path','metadata'):
        if spec['completeness'].get(key):spec['completeness'][key]=str((root/spec['completeness'][key]).resolve())
    nodes,weights=normal_nodes(constraint['slope_mean']+offset,constraint['slope_sigma'],order)
    configs=[];priors={}
    for i,(slope,weight) in enumerate(zip(nodes,weights)):
        description=f'{cluster}: explicitly assumed indirect global-MF estimator model; slope={slope:g}, delta={offset:g}; not a measured escaping MF or validated absolute mass-loss rate.'
        pop=measured_population(isochrone,slope,age=imeta['age_years'],mh=imeta['metallicity_mh'],source=description)
        pop.metadata.update(global_mf_constraint=constraint,mass_function_status='indirect_estimator_model_assumption_not_measured',
            escape_slope_offset=offset,quadrature_prior_weight=float(weight),physical_rate_measurement_validated=False,
            population_policy_source_sha256=sha256_file(Path(__file__)))
        pop.save(out/f'population_{i}')
        scenario=dict(empirical_mf_slope=float(slope),escape_slope_offset=float(offset),
            parallax_zero_point_mas=spec['parallax_zero_point_mas'],mf_constraint_kind=constraint['constraint_kind'])
        s=dict(spec,population=str((out/f'population_{i}').resolve()),scenario=scenario,source_description=description)
        s.pop('lookback_window_myr',None)
        p=out/f'response_{i}.json';p.write_text(json.dumps(s,indent=2),encoding='utf-8');configs.append(str(p.resolve()))
        priors[json.dumps(scenario,sort_keys=True)]=float(weight)
    report=dict(cluster=cluster,constraint=constraint,escape_slope_offset=offset,quadrature_order=order,
        scenario_weights=priors,response_configs=configs,source_sha256=sha256_file(Path(__file__)),
        reference_config_sha256=sha256_file(reference_config),isochrone_sha256=sha256_file(isochrone),
        interpretation='Conditional indirect-estimator model; direct MF measurement and a validated absolute-rate interpretation are unavailable.')
    (out/'population_prior.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report
