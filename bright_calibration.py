"""Extend a prepared calibration using an independently queried bright census."""
import argparse,json,os,shutil,subprocess,sys
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from provenance import sha256_file


def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,d):Path(p).write_text(json.dumps(d,indent=2),encoding='utf-8')


def prepare(repair,cluster,bright,out):
    repair,bright,out=map(lambda p:Path(p).resolve(),(repair,bright,out));info=read(repair/'repair_manifest.json')
    for path,sha in info['frozen_sha256'].items():
        if sha256_file(path)!=sha:raise ValueError('Frozen source changed: '+path)
    if out.exists():raise FileExistsError('Choose a fresh support-repair version')
    source=Path(info['source'])/cluster;selection=read(source/'selection/selection.json')
    faint=source/'calibration_inputs/parent_random10.csv.gz'
    if not faint.exists():
        candidate=Path(__file__).resolve().parent/'outputs/nine_cluster_real_inputs_v1'/cluster/'parent_random10.csv.gz'
        expected=read(source/'science_calibration/bundle/calibration_manifest.json')['parent_source']['original_parent_sha256']
        if not candidate.exists() or sha256_file(candidate)!=expected:raise ValueError('Original uncut faint parent unavailable or checksum mismatch')
        faint=candidate
    fm=read(faint.with_suffix('.json'));bm=read(bright.with_suffix('.json'))
    ra,dec=selection['center_ra_deg'],selection['center_dec_deg'];radius=selection['field_radius_deg']+.05
    from selection_model import parent_query
    if bm['query_adql']!=parent_query(ra,dec,radius,0.,12.,1):raise ValueError('Bright census query differs from the required geometry and G strata')
    if not (fm['g_min']==12 and fm['g_max']==21 and fm['random_sampling_modulus']==10 and fm['radius_deg']==radius and np.allclose(fm['center'],[ra,dec])):
        raise ValueError('Unsupported faint sampling design')
    for p,m in ((faint,fm),(bright,bm)):
        if m['sha256']!=sha256_file(p) or m['query_status']!='OK' or not m['parent_is_unfiltered']:raise ValueError('Unverified component parent')
    out.mkdir(parents=True);dest=out/cluster;inp=dest/'calibration_inputs';inp.mkdir(parents=True)
    for suffix in ('dat','json'):shutil.copy2(source/'calibration_inputs'/('isochrone.'+suffix),inp/('isochrone.'+suffix))
    a,b=pd.read_csv(faint,low_memory=False),pd.read_csv(bright,low_memory=False)
    if not ((a.phot_g_mean_mag>=12)&(a.phot_g_mean_mag<21)).all() or not ((b.phot_g_mean_mag>=0)&(b.phot_g_mean_mag<12)).all():raise ValueError('Component magnitude strata overlap or exceed their declared support')
    combined=pd.concat([a,b],ignore_index=True)
    if combined.source_id.duplicated().any():raise ValueError('Repeated Gaia source across disjoint strata')
    parent=inp/'parent_stratified.csv.gz';combined.to_csv(parent,index=False,compression='gzip')
    query=parent_query(ra,dec,radius,0.,21.,1)+ ' AND ((phot_g_mean_mag < 12) OR (phot_g_mean_mag >= 12 AND MOD(random_index,10) = 0))'
    metadata=dict(parent_is_unfiltered=True,synthetic_calibration=False,query_status='OK',query_adql=query,
        source_description='Gaia DR3 sky/G-only stratified calibration parent: full census below G=12, independent 10% sample at G=12..21; no quality or membership cuts. Sampling changes only at a calibration-bin edge.',
        n_sources=len(combined),sha256=sha256_file(parent),radius_deg=radius,center=[ra,dec],g_min=0.,g_max=21.,
        sampling_strata=[dict(g_range=[0,12],fraction=1.),dict(g_range=[12,21],fraction=.1)],
        components=[dict(metadata=fm,original_path=str(faint)),dict(metadata=bm,original_path=str(bright))],
        uncertainty='Within each nonoverlapping G stratum sampling is independent of availability and quality flags; raw conditional binomial counts are used. No global sampling fraction is assigned.')
    dump(parent.with_suffix('.json'),metadata)
    del a,b,combined
    code=out/'frozen_code';shutil.copytree(info['snapshot'],code,ignore=shutil.ignore_patterns('__pycache__','outputs'))
    shutil.copy2(Path(__file__),code/'bright_calibration.py')
    p=code/'prepare_science_calibration.py';text=p.read_text(encoding='utf-8');needle='np.array([12.,18.,19.,20.,20.5,21.])';assert text.count(needle)==1
    p.write_text(text.replace(needle,'np.array([0.,12.,18.,19.,20.,20.5,21.])'),encoding='utf-8')
    os.environ['GC_TIDAL_STUDY']=info['study']
    master=pd.read_csv(Path(info['study'])/'catalogs/master_clusters.csv').set_index('id');bundle=dest/'science_calibration/bundle'
    command=[sys.executable,str(code/'prepare_science_calibration.py'),'--isochrone',str(inp/'isochrone.dat'),'--parent',str(parent),
        '--selection',str(source/'selection'),'--study',info['study'],'--ra',str(ra),'--dec',str(dec),'--radius',str(radius),
        '--ebv',str(master.loc[cluster,'ebv']),'--out',str(bundle)]
    with (out/'calibration.log').open('w',encoding='utf-8') as log:subprocess.run(command,cwd=code,stdout=log,stderr=subprocess.STDOUT,check=True)
    # Subsequent imports must come from the new snapshot, not from root modules.
    config=dict(cluster=cluster,source=info['source'],snapshot=str(code),study=info['study'],parent_repair=str(repair),
        calibration_inputs_root=str(inp),reference=str(bundle/'tilt_plus0_zp_minus0.017.json'),
        created_utc=datetime.now(timezone.utc).isoformat(),physical_rate_measurement_validated=False)
    dump(out/'bright_repair.json',config)
    print(json.dumps(config,indent=2))


def finalize(out):
    out=Path(out).resolve();settings=read(out/'bright_repair.json');code=Path(settings['snapshot']);cluster=settings['cluster']
    os.environ['GC_TIDAL_STUDY']=settings['study'];sys.path.insert(0,str(code))
    from empirical_population import prepare as population
    from population_model import Population
    from selection_model import QualityMap
    from response_group import run_group
    from validate_science import command,check_response
    source=Path(settings['source'])/cluster;inp=Path(settings['calibration_inputs_root']);dest=out/cluster
    quality=QualityMap.load(Path(settings['reference']).parent/'quality')
    # This extension does not coarsen cells or lower their support threshold.
    if min(quality.arrays['n_parent'].min(),quality.arrays['n_available'].min())<quality.minimum:
        raise ValueError('Bright extension still has sparse cells; review calibration before running')
    # The original MF table is copied from the versioned parent repair supplied
    # during preparation, never downloaded or silently replaced at finalization.
    original=Path(settings['parent_repair'])/'mf_source.html'
    shutil.copy2(original,out/'mf_source.html')
    reports=[]
    for offset in (0.,-.5,-1.):
        folder=dest/'empirical_repair_v1'/('offset_'+f'{offset:g}'.replace('-','minus').replace('.','p'))
        prior=population(cluster,inp/'isochrone.dat',settings['reference'],out/'mf_source.html',folder/'population_bundle',offset=offset,order=5)
        for response in prior['response_configs']:
            pop=Population.load(read(response)['population'])
            if not (pop.table.number_per_msun>0).all() or not np.isclose(np.sum(pop.table.number_per_msun*pop.table.current_mass_msun),1.,rtol=1e-12):raise ValueError('Invalid population normalization')
        reports.append(dict(offset=offset,prior=str(folder/'population_bundle/population_prior.json'),nodes=5))
    c=read(source/'coupled_dynamics.json');c.update(duration_myr=2.,n_particles=16,release_sampling='stratified',python_integrator='cartesian_batch',batch_size=4,integration_rtol=1e-11,integration_atol=1e-13)
    smoke=out/'preflight';smoke.mkdir();dump(smoke/'config.json',c)
    command('stream_model.py',['--config',smoke/'config.json','--out',smoke/'simulation'],smoke/'integration.log')
    middle=read(reports[0]['prior'])['response_configs'][2]
    run_group(source/'selection',smoke/'simulation',[middle],[smoke/'response'],study=Path(settings['study']))
    response=check_response(smoke/'response')
    preflight=dict(stage='passed_small_real_input_response',tracers=16,duration_myr=2.,
        min_parent=int(quality.arrays['n_parent'].min()),min_available=int(quality.arrays['n_available'].min()),
        response_template_sha256=response['template_sha256'],physical_rate_measurement_validated=False,
        scope='Small real-input software/observation check; full-duration spatial, donor, template convergence and scientific recovery checks remain necessary.')
    dump(out/'preflight.json',preflight)
    state=dict(stage='ready_for_conditional_inference',reference=settings['reference'],population_checks=reports,
        calibration_inputs_root=str(inp),representative_response_passed=True,preflight=preflight,
        updated_utc=datetime.now(timezone.utc).isoformat(),physical_rate_measurement_validated=False)
    dump(out/'input_status.json',{cluster:state})
    files=[*code.glob('*.py'),*list((code/'configs').rglob('*.json')),*[p for p in dest.rglob('*') if p.is_file()],out/'mf_source.html',
        source/'coupled_dynamics.json',source/'selection/selection.json',source/'selection/maps.npz',Path(settings['study'])/'step3_tail_search.py']
    manifest=dict(schema_version=1,source=settings['source'],snapshot=str(code),study=settings['study'],clusters=[cluster],
        created_utc=datetime.now(timezone.utc).isoformat(),physical_rate_measurement_validated=False,
        policy='bright_census_plus_random_faint_parent_v1; original spatial/colour/G>=12 bins and count threshold retained',
        frozen_sha256={str(p.resolve()):sha256_file(p) for p in files})
    dump(out/'repair_manifest.json',manifest);print(json.dumps(preflight,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--finalize',action='store_true')
    for key in ('repair','bright'):p.add_argument('--'+key,type=Path)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--cluster');a=p.parse_args()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    if a.finalize:finalize(a.out)
    else:
        if a.repair is None or a.bright is None or a.cluster is None:p.error('--repair, --bright and --cluster required')
        prepare(a.repair,a.cluster,a.bright,a.out)
