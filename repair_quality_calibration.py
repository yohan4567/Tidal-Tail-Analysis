"""Prepare and validate a versioned sparse-quality repair without changing live jobs.

An isolated source snapshot enables pooled quality maps. Existing production
modules, calibrations, dynamics and posteriors are never rewritten.
"""
import argparse
from datetime import datetime, timezone
import json,os,shutil,subprocess,sys
from pathlib import Path
import pandas as pd
from provenance import sha256_file
from quality_repair import validate_pooling
from selection_model import QualityMap
from selection_bridge import DEFAULT_STUDY
from validate_science import absolute_spec

HERE=Path(__file__).resolve().parent


def dump(path,data):
    Path(path).write_text(json.dumps(data,indent=2),encoding='utf-8')


def verify_snapshot(repair):
    info=json.loads((repair/'repair.json').read_text(encoding='utf-8'))
    for path,digest in info['frozen_sha256'].items():
        if sha256_file(path)!=digest:raise ValueError('Frozen repair input changed: '+path)
    return info


def prepare(root,cluster,out):
    root,out=Path(root).resolve(),Path(out).resolve();dest=root/cluster
    if out.exists():raise FileExistsError('Use a fresh repair directory')
    bundle=dest/'science_calibration/bundle';quality=QualityMap.load(bundle/'quality')
    holdout=bundle/'quality_holdout_cells.csv'
    candidate,validation=validate_pooling(quality,pd.read_csv(holdout))
    if not validation['screening_passed']:raise ValueError('Pooled quality failed held-out field screening; obtain more parent data')
    source_grid=dest/'empirical_direct_v1/offset_0/grid'
    plan=json.loads((source_grid/'plan.json').read_text())
    for rate in plan['rates']:
        for seed in plan['seeds']:
            sim=source_grid/f'rate_{rate:g}_seed_{seed}'/'simulation'
            if not (sim/'run.json').exists():raise ValueError('Repair requires completed dynamics; missing '+str(sim))
    out.mkdir(parents=True);candidate.save(out/'quality')
    dump(out/'holdout_validation.json',validation)
    spec=absolute_spec(bundle/'tilt_plus0_zp_minus0.017.json')
    spec.update(quality=str(out/'quality'),cache_dir=str(out/'response_cache'),
        quality_repair=dict(policy=candidate.metadata['quality_pooling']['policy'],
            assumption=candidate.metadata['quality_pooling']['assumption']))
    dump(out/'reference_response.json',spec)
    code=out/'frozen_code';code.mkdir()
    for p in HERE.glob('*.py'):shutil.copy2(p,code/p.name)
    shutil.copytree(HERE/'configs',code/'configs')
    if (HERE/'.galpyrc').exists():shutil.copy2(HERE/'.galpyrc',code/'.galpyrc')
    # These explicit adaptations are confined to the new source snapshot.
    path=code/'selection_model.py';text=path.read_text(encoding='utf-8')
    needle='        return cls(arrays,meta,draw_seed)'
    assert text.count(needle)==1
    text=text.replace(needle,"        if meta.get('quality_pooling'):\n            from quality_repair import SpatialPooledQualityMap\n            return SpatialPooledQualityMap(arrays,meta,draw_seed)\n"+needle)
    path.write_text(text,encoding='utf-8')
    path=code/'response_model.py';text=path.read_text(encoding='utf-8')
    needle="'population_model.py','joint_grid.py')";assert text.count(needle)==1
    path.write_text(text.replace(needle,"'population_model.py','joint_grid.py','quality_repair.py')"),encoding='utf-8')
    path=code/'run_joint_study.py';text=path.read_text(encoding='utf-8')
    needle="'selection_model.py','validate_window_reuse.py'";assert text.count(needle)==1
    path.write_text(text.replace(needle,"'selection_model.py','quality_repair.py','validate_window_reuse.py'"),encoding='utf-8')
    baseline=dest/'empirical_direct_v1/offset_0/population_bundle/response_0.json'
    info=dict(schema_version=1,created_utc=datetime.now(timezone.utc).isoformat(),cluster=cluster,
        source_root=str(root),source_grid=str(source_grid),selection=str(dest/'selection'),
        config=str(dest/'coupled_dynamics.json'),isochrone=str(dest/'calibration_inputs/isochrone.dat'),
        html=str(HERE/'outputs/empirical_mf_source/parameter.html'),reference=str(out/'reference_response.json'),
        study=str(DEFAULT_STUDY.resolve()),
        validation_population_config=str(baseline),source_snapshot=str(code),
        holdout_screening=validation['screening_passed'],pooling=candidate.metadata['quality_pooling'],
        physical_rate_measurement_validated=False,
        scope='Conditional spatial pooling; no new reference stars or reduced min_parent; original dynamics reused only after signature checks.')
    files=[*code.glob('*.py'),*list((code/'configs').rglob('*.json')),out/'quality/quality.json',out/'quality/quality.npz',
        out/'reference_response.json',out/'holdout_validation.json',Path(info['config']),
        Path(info['selection'])/'selection.json',Path(info['selection'])/'maps.npz',Path(info['isochrone']),
        Path(info['isochrone']).with_suffix('.json'),Path(info['html']),bundle/'quality/quality.json',bundle/'quality/quality.npz',holdout]
    info['frozen_sha256']={str(p.resolve()):sha256_file(p) for p in files}
    dump(out/'repair.json',info)
    print(json.dumps(dict(repair=str(out),snapshot=str(code),pooled_strata=len(info['pooling']['pooled_strata']),
        minimum_pooled_n=int(candidate.group_n.min()),holdout_max_abs_z=validation['max_abs_predictive_z']),indent=2))


def validate(repair):
    repair=Path(repair).resolve();info=verify_snapshot(repair)
    if HERE!=Path(info['source_snapshot']):raise ValueError('Run validation from the frozen_code script recorded in repair.json')
    os.environ['GC_TIDAL_STUDY']=info['study']
    from response_group import run_group
    from validate_science import check_response
    spec=absolute_spec(info['validation_population_config'])
    spec.update(quality=str(repair/'quality'),cache_dir=str(repair/'validation_cache'),quality_repair=info['pooling'])
    cfg=repair/'validation_response.json';dump(cfg,spec)
    plan=json.loads((Path(info['source_grid'])/'plan.json').read_text())
    simulation=Path(info['source_grid'])/f"rate_{plan['rates'][0]:g}_seed_{plan['seeds'][0]}"/'simulation'
    out=repair/'validation_response'
    state=repair/'validation_status.json';dump(state,dict(stage='running',simulation=str(simulation),physical_rate_measurement_validated=False))
    try:
        run_group(Path(info['selection']),simulation,[cfg],[out],study=Path(info['study']))
        check_response(out)
        result=json.loads((out/'response.json').read_text())
        dump(state,dict(stage='completed_representative_response',simulation=str(simulation),
            expected_stars_in_fit_region=result['expected_stars_in_fit_region'],
            observation_draws=spec['observation_draws'],calibration_draws=spec['calibration_draws'],
            physical_rate_measurement_validated=False))
    except Exception as exc:
        dump(state,dict(stage='failed',error=str(exc),physical_rate_measurement_validated=False));raise


def run(repair,workers):
    repair=Path(repair).resolve();info=verify_snapshot(repair)
    if HERE!=Path(info['source_snapshot']):raise ValueError('Run the frozen_code script recorded in repair.json')
    os.environ['GC_TIDAL_STUDY']=info['study']
    validation=json.loads((repair/'validation_status.json').read_text())
    if validation['stage']!='completed_representative_response':raise ValueError('Representative real-input response must pass first')
    from empirical_population import prepare as population
    from run_joint_study import run as grid
    from bayesian_grid import run as infer
    plan=json.loads((Path(info['source_grid'])/'plan.json').read_text());reuse=Path(info['source_grid']);reports=[]
    state=repair/'run_status.json'
    for offset in (0.,-.5,-1.):
        folder=repair/('offset_'+f'{offset:g}'.replace('-','minus').replace('.','p'))
        bundle=folder/'population_bundle';prior_path=bundle/'population_prior.json'
        if prior_path.exists():prior=json.loads(prior_path.read_text())
        else:prior=population(info['cluster'],info['isochrone'],info['reference'],info['html'],bundle,offset=offset,order=5)
        dump(state,dict(stage='running_conditional_grid',offset=offset,workers=workers,physical_rate_measurement_validated=False))
        grid(info['config'],prior['response_configs'],info['selection'],folder/'grid',rates=plan['rates'],times=plan['times_myr'],
            seeds=plan['seeds'],interval_myr=plan['release_interval_myr'],workers=workers,response_workers=workers,
            cluster_age_myr=plan['cluster_age_myr'],pixel_factor=plan.get('pixel_factor',1),
            integrator=plan.get('integrator','cartesian'),batch_size=plan.get('batch_size',32),reuse_grid=reuse)
        infer(folder/'grid',folder/'posterior',prior_path);reports.append(str(folder/'posterior/posterior_summary.json'))
        if offset==0:reuse=folder/'grid'
    dump(state,dict(stage='completed_conditional_posteriors_validation_pending',reports=reports,
        pooling_assumption=info['pooling']['assumption'],physical_rate_measurement_validated=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare');a.add_argument('--root',type=Path,required=True);a.add_argument('--cluster',required=True);a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('validate');a.add_argument('--repair',type=Path,required=True)
    a=sub.add_parser('run');a.add_argument('--repair',type=Path,required=True);a.add_argument('--workers',type=int,default=12)
    args=p.parse_args()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    if args.action=='prepare':prepare(args.root,args.cluster,args.out)
    elif args.action=='validate':validate(args.repair)
    else:
        if not 1<=args.workers<=12:p.error('Use 1–12 workers')
        run(args.repair,args.workers)
