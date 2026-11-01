"""Run a source-frozen repair grid and likelihood comparison, without R/T posterior."""
import argparse,json,os,sys
from pathlib import Path
from datetime import datetime,timezone
from provenance import sha256_file


def run(job,out,workers=12,validate_only=False):
    job,out=Path(job),Path(out);j=json.loads(job.read_text());code=Path(j['code']);out.mkdir(parents=True,exist_ok=True)
    for path,sha in j['frozen_sha256'].items():
        if sha256_file(path)!=sha:raise ValueError('Frozen priority input changed: '+path)
    os.environ['GC_TIDAL_STUDY']=j['study'];sys.path.insert(0,str(code))
    from empirical_population import prepare
    from run_joint_study import run as grid
    from production_study import RATES,TIMES,SEEDS
    from likelihood_comparison import compare
    state=out/'run_status.json';reuse=Path(j['reuse_grid']) if j.get('reuse_grid') else None;reports=[]
    priors={}
    from population_model import Population
    import numpy as np
    support=None
    for offset in (0.,-.5,-1.):
        folder=out/('offset_'+f'{offset:g}'.replace('-','minus').replace('.','p'));bundle=folder/'population_bundle';prior_path=bundle/'population_prior.json'
        prior=json.loads(prior_path.read_text()) if prior_path.exists() else prepare(j['cluster'],j['isochrone'],j['reference'],j['html'],bundle,offset=offset,order=5)
        for cfg in prior['response_configs']:
            pop=Population.load(json.loads(Path(cfg).read_text())['population']);t=pop.table
            if not (t.number_per_msun>0).all() or not np.isclose(np.sum(t.number_per_msun*t.current_mass_msun),1.,rtol=1e-12):raise ValueError('Invalid population support or mass normalization')
            values=t.drop(columns='number_per_msun')
            if support is None:support=values
            elif not support.equals(values):raise ValueError('Photometric support differs across population scenarios')
        priors[offset]=(prior_path,prior)
    validation=out/'preflight.json'
    if not validation.exists():
        from validate_science import command,check_response
        from response_group import run_group
        simulation=None
        if reuse is not None:
            old=json.loads((reuse/'plan.json').read_text());simulation=reuse/f"rate_{old['rates'][0]:g}_seed_{old['seeds'][0]}"/'simulation'
            if not (simulation/'run.json').exists():simulation=None
        if simulation is None:
            c=json.loads(Path(j['config']).read_text());c.update(duration_myr=2.,n_particles=16,release_sampling='stratified',python_integrator='cartesian_batch',batch_size=4,integration_rtol=1e-11,integration_atol=1e-13)
            cfg=out/'preflight_config.json';cfg.write_text(json.dumps(c,indent=2));simulation=out/'preflight_simulation'
            command('stream_model.py',['--config',cfg,'--out',simulation],out/'preflight_integration.log')
        response=out/'preflight_response';run_group(Path(j['selection']),simulation,[priors[0.][1]['response_configs'][2]],[response],study=Path(j['study']));check_response(response)
        validation.write_text(json.dumps(dict(stage='passed_representative_response',job_sha256=sha256_file(job),simulation=str(simulation),population_nodes_checked=15,full_physical_validation=False),indent=2))
    elif json.loads(validation.read_text())['job_sha256']!=sha256_file(job):raise ValueError('Preflight job changed')
    if validate_only:return
    for offset in (0.,-.5,-1.):
        folder=out/('offset_'+f'{offset:g}'.replace('-','minus').replace('.','p'));bundle=folder/'population_bundle';prior_path=bundle/'population_prior.json'
        prior_path,prior=priors[offset]
        state.write_text(json.dumps(dict(stage='running_model_grid',cluster=j['cluster'],offset=offset,workers=workers,updated_utc=datetime.now(timezone.utc).isoformat(),target_parameter_posterior_computed=False),indent=2))
        grid(j['config'],prior['response_configs'],j['selection'],folder/'grid',rates=RATES,times=TIMES,seeds=SEEDS,
            interval_myr=.25,workers=workers,response_workers=workers,cluster_age_myr=j['cluster_age_myr'],pixel_factor=4,
            integrator='cartesian_batch',batch_size=128,reuse_grid=reuse)
        state.write_text(json.dumps(dict(stage='integrating_nuisance_parameters',cluster=j['cluster'],offset=offset,updated_utc=datetime.now(timezone.utc).isoformat()),indent=2))
        compare(folder/'grid',folder/'comparison',prior_path);reports.append(str(folder/'comparison/comparison.json'))
        if offset==0:reuse=folder/'grid'
    state.write_text(json.dumps(dict(stage='completed_model_comparisons_validation_pending',cluster=j['cluster'],reports=reports,
        target_parameter_posterior_computed=False,physical_rate_measurement_validated=False,updated_utc=datetime.now(timezone.utc).isoformat()),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--job',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--workers',type=int,default=12);p.add_argument('--validate-only',action='store_true');a=p.parse_args()
    if not 1<=a.workers<=12:p.error('Use 1–12 workers on the local PC')
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    run(a.job,a.out,a.workers,a.validate_only)
