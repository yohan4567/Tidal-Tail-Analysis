"""Isolated, auditable repair of seven survey inputs; never rewrite live results."""
import argparse,json,os,shutil,sys
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from provenance import sha256_file

HERE=Path(__file__).resolve().parent
CLUSTERS=('NGC3201','NGC362','NGC6121','NGC6171','NGC6397','NGC6266','NGC7089')


def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def dump(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2),encoding='utf-8')


def verify(out):
    info=read(out/'repair_manifest.json')
    for file,digest in info['frozen_sha256'].items():
        if sha256_file(file)!=digest:raise ValueError('Frozen repair input changed: '+file)
    if HERE!=Path(info['snapshot']):raise ValueError('Execute the frozen snapshot script')
    os.environ['GC_TIDAL_STUDY']=info['study']
    return info


def bootstrap(source,out):
    from selection_bridge import DEFAULT_STUDY
    source,out=Path(source).resolve(),Path(out).resolve()
    if out.exists():raise FileExistsError('Choose a new versioned repair directory')
    out.mkdir(parents=True);code=out/'frozen_code';code.mkdir()
    for p in HERE.glob('*.py'):shutil.copy2(p,code/p.name)
    shutil.copytree(HERE/'configs',code/'configs')
    if (HERE/'.galpyrc').exists():shutil.copy2(HERE/'.galpyrc',code/'.galpyrc')
    # Numerical integration sources remain byte-for-byte identical, preserving
    # authentic checkpoint identities. Route the process through an I/O wrapper.
    p=code/'validate_science.py';text=p.read_text(encoding='utf-8')
    needle='    log=Path(log); log.parent.mkdir(parents=True,exist_ok=True)'
    assert text.count(needle)==1
    p.write_text(text.replace(needle,"    if script=='stream_model.py':script='checkpoint_io.py'\n"+needle),encoding='utf-8')
    p=code/'empirical_population.py';text=p.read_text(encoding='utf-8')
    start=text.index('    exponent=slope+1.;ratio=np.log(proxy[1:]/proxy[:-1])')
    stop=text.index("    t['number_per_msun']",start)
    text=text[:start]+"    from population_quadrature_repair import mass_quadrature\n    weights,quadrature=mass_quadrature(proxy,slope,t.number_per_msun.to_numpy())\n"+text[stop:]
    needle="meta=dict(base.metadata,weight_mode='empirical_present_day_ms_power_law_v1'"
    assert text.count(needle)==1
    text=text.replace(needle,"meta=dict(base.metadata,mass_quadrature=quadrature,weight_mode='empirical_present_day_ms_power_law_v2'")
    p.write_text(text,encoding='utf-8')
    p=code/'run_joint_study.py';text=p.read_text(encoding='utf-8')
    needle="'selection_model.py','validate_window_reuse.py'";assert text.count(needle)==1
    p.write_text(text.replace(needle,"'selection_model.py','population_quadrature_repair.py','checkpoint_io.py','validate_science.py','validate_window_reuse.py'"),encoding='utf-8')
    # Run the adaptive downloader through the existing survey CLI, avoiding a
    # circular module import or any change to its tiling query construction.
    p=code/'run_remaining_survey.py';text=p.read_text(encoding='utf-8')
    p.write_text(text.replace("execute('fetch_parent_tiles.py'","execute('adaptive_parent.py'"),encoding='utf-8')
    p=code/'prepare_science_calibration.py';text=p.read_text(encoding='utf-8')
    old="str((Path(__file__).resolve().parent/'outputs/cache/gaia/allsky_uniq_10.fits').resolve())"
    assert text.count(old)==1
    p.write_text(text.replace(old,repr(str((HERE/'outputs/cache/gaia/allsky_uniq_10.fits').resolve()))),encoding='utf-8')
    initial=[]
    for cid in CLUSTERS:
        dest=out/cid;dest.mkdir();old=source/cid
        for name in ('selection','calibration_inputs'):
            if name=='selection':shutil.copytree(old/name,dest/name)
            else:
                (dest/name).mkdir()
                for suffix in ('dat','json'):shutil.copy2(old/name/('isochrone.'+suffix),dest/name/('isochrone.'+suffix))
        for name in ('adopted_dynamics.json','coupled_dynamics.json'):
            if (old/name).exists():shutil.copy2(old/name,dest/name)
        initial.extend(p for p in dest.rglob('*') if p.is_file())
    html=HERE/'outputs/empirical_mf_source/parameter.html'
    shutil.copy2(html,out/'mf_source.html');initial.append(out/'mf_source.html')
    files=[*code.glob('*.py'),*list((code/'configs').rglob('*.json')),*initial,
        HERE/'outputs/cache/gaia/allsky_uniq_10.fits',DEFAULT_STUDY/'step3_tail_search.py']
    info=dict(schema_version=1,created_utc=datetime.now(timezone.utc).isoformat(),source=str(source),snapshot=str(code),
        study=str(DEFAULT_STUDY.resolve()),clusters=list(CLUSTERS),
        policy='input_repairs_v1; tied mass quadrature + overflow-safe parent + bounded checkpoint replace retries',
        physical_rate_measurement_validated=False,frozen_sha256={str(p.resolve()):sha256_file(p) for p in files})
    dump(out/'repair_manifest.json',info)
    print(json.dumps(dict(repair=str(out),snapshot=str(code)),indent=2))


def population_check(out,cid,reference,info):
    from empirical_population import prepare
    from population_model import Population
    dest=out/cid;reports=[]
    for offset in (0.,-.5,-1.):
        folder=dest/'empirical_repair_v1'/('offset_'+f'{offset:g}'.replace('-','minus').replace('.','p'))
        prior_path=folder/'population_bundle/population_prior.json'
        prior=read(prior_path) if prior_path.exists() else prepare(cid,dest/'calibration_inputs/isochrone.dat',reference,
            out/'mf_source.html',prior_path.parent,offset=offset,order=5)
        for spec in prior['response_configs']:
            pop=Population.load(read(spec)['population']);t=pop.table
            if not np.isfinite(t.number_per_msun).all() or not (t.number_per_msun>0).all():raise ValueError('Population support invalid')
            if not np.isclose(np.sum(t.current_mass_msun*t.number_per_msun),1.,rtol=1e-12):raise ValueError('Population mass normalization failed')
        reports.append(dict(offset=offset,prior=str(prior_path),nodes=len(prior['response_configs'])))
    return reports


def prepare_inputs(out):
    info=verify(out);source=Path(info['source']);status=read(out/'input_status.json') if (out/'input_status.json').exists() else {}
    os.environ['GC_PARENT_REUSE_DIRS']=json.dumps([str(source/c/'calibration_inputs/parent_random10.csv.gz.parts') for c in CLUSTERS])
    from run_remaining_survey import run as prepare
    for cid in CLUSTERS:
        if status.get(cid,{}).get('stage')=='ready_for_conditional_inference':continue
        state=dict(stage='checking_inputs',updated_utc=datetime.now(timezone.utc).isoformat(),physical_rate_measurement_validated=False)
        status[cid]=state;dump(out/'input_status.json',status)
        try:
            reference=source/cid/'science_calibration/bundle/tilt_plus0_zp_minus0.017.json'
            if cid not in ('NGC362','NGC3201'):
                state['stage']='downloading_or_calibrating';dump(out/'input_status.json',status)
                prepare(out,HERE/'configs/survey_populations_v1.json',workers=1,response_workers=1,prepare_workers=1,clusters=[cid],prepare_only=True)
                job=read(out/cid/'remaining_job.json')
                if job['stage']!='ready_conditional_grid':raise RuntimeError(job.get('error',job.get('reason',job['stage'])))
                reference=out/cid/'science_calibration/bundle/tilt_plus0_zp_minus0.017.json'
            reports=population_check(out,cid,reference,info)
            state.update(stage='ready_for_conditional_inference',reference=str(reference),population_checks=reports,
                interpretation='Positive support and unit represented mass checked for all 15 MF nodes; full response support and physical-rate validation still required.')
        except Exception as exc:state.update(stage='failed',error=str(exc))
        state['updated_utc']=datetime.now(timezone.utc).isoformat();dump(out/'input_status.json',status)
        print(cid,state['stage'],state.get('error',''),flush=True)


def run_repaired(out,workers):
    info=verify(out);source=Path(info['source']);status=read(out/'input_status.json')
    from run_joint_study import run as grid
    from bayesian_grid import run as infer
    from production_study import RATES,TIMES,SEEDS
    rules=read(HERE/'configs/survey_populations_v1.json')['clusters']
    for cid,extra in read(HERE/'configs/population_supplements_v1.json')['clusters'].items():rules[cid]=dict(rules.get(cid,{}),**extra)
    for cid in CLUSTERS:
        if status.get(cid,{}).get('stage')!='ready_for_conditional_inference':continue
        state_path=out/cid/'inference_status.json';old=read(state_path) if state_path.exists() else {}
        if old.get('stage')=='completed_conditional_posteriors_validation_pending':continue
        state=dict(stage='running',workers=workers,physical_rate_measurement_validated=False)
        dump(state_path,state)
        try:
            reuse=source/cid/'empirical_direct_v1/offset_0/grid'
            if not (reuse/'plan.json').exists():reuse=None
            for item in status[cid]['population_checks']:
                prior_path=Path(item['prior']);prior=read(prior_path);folder=prior_path.parent.parent
                state.update(offset=item['offset'],updated_utc=datetime.now(timezone.utc).isoformat());dump(state_path,state)
                grid(out/cid/'coupled_dynamics.json',prior['response_configs'],out/cid/'selection',folder/'grid',
                    rates=RATES,times=TIMES,seeds=SEEDS,interval_myr=.25,workers=workers,response_workers=workers,
                    cluster_age_myr=rules[cid]['age_gyr']*1000,pixel_factor=4,integrator='cartesian_batch',batch_size=128,reuse_grid=reuse)
                infer(folder/'grid',folder/'posterior',prior_path)
                if item['offset']==0:reuse=folder/'grid'
            state['stage']='completed_conditional_posteriors_validation_pending'
        except Exception as exc:state.update(stage='failed',error=str(exc))
        state['updated_utc']=datetime.now(timezone.utc).isoformat();dump(state_path,state)
        print(cid,state['stage'],state.get('error',''),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('bootstrap','prepare','run'))
    p.add_argument('--out',type=Path,required=True);p.add_argument('--source',type=Path)
    p.add_argument('--workers',type=int,default=12);a=p.parse_args()
    if not 1<=a.workers<=12:p.error('Use 1–12 workers')
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    if a.action=='bootstrap':
        if a.source is None:p.error('--source required for bootstrap')
        bootstrap(a.source,a.out)
    elif a.action=='prepare':prepare_inputs(a.out.resolve())
    else:run_repaired(a.out.resolve(),a.workers)
