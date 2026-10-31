"""Run empirical-MF response grids and marginalized posteriors on verified dynamics.

Separate escape-transfer offsets are sensitivity experiments, not equally weighted
population models. Existing frozen grids remain unchanged.
"""
import argparse,json,os
from datetime import datetime,timezone
from pathlib import Path
from empirical_population import prepare
from run_joint_study import run as run_grid
from bayesian_grid import run as infer
from provenance import sha256_file

HERE=Path(__file__).resolve().parent


def run(root,clusters=('Pal5','NGC288'),offsets=(0.,-.5,-1.),order=5,response_workers=2,study_label='empirical_mf_v1'):
    root=Path(root).resolve();html=HERE/'outputs/empirical_mf_source/parameter.html'
    if Path(study_label).name!=study_label or study_label in ('.','..'):raise ValueError('Study label must be a directory basename')
    if not html.exists():
        from network import open_url
        html.parent.mkdir(parents=True,exist_ok=True)
        with open_url('https://people.smp.uq.edu.au/HolgerBaumgardt/globular/parameter.html') as response:
            content=response.read()
        with html.open('xb') as target:target.write(content)
    for cid in clusters:
        old=root/cid/'production_grid_shared_v1'
        if not (old/'scores/joint_grid.json').exists():raise ValueError('Need a completed verified source grid: '+cid)
    for cid in clusters:
        old=root/cid/'production_grid_shared_v1';plan=json.loads((old/'plan.json').read_text())
        selection=root/cid/'selection';state_path=root/cid/'empirical_job.json'
        def update(stage,**kw):
            data=dict(cluster=cid,stage=stage,updated_utc=datetime.now(timezone.utc).isoformat(),
                physical_rate_measurement_validated=False,**kw)
            temp=state_path.with_suffix('.tmp');temp.write_text(json.dumps(data,indent=2),encoding='utf-8');temp.replace(state_path)
            print(cid,stage,kw.get('offset',''),flush=True)
        iso=HERE/'outputs/Pal5_science_calibration/isochrone.dat' if cid=='Pal5' else root/cid/'calibration_inputs/isochrone.dat'
        cfg=old/'configs'/f'rate_{plan["rates"][0]:g}_seed_{plan["seeds"][0]}.json'
        reference=old/'configs'/f'rate_{plan["rates"][0]:g}_seed_{plan["seeds"][0]}_pop_0_T_{max(plan["times_myr"]):g}.json'
        reports=[]
        for offset in offsets:
            label=f'offset_{offset:g}'.replace('-','minus').replace('.','p')
            dest=root/cid/study_label/label;bundle=dest/'population_bundle';grid=dest/'grid'
            update('preparing_measured_population',offset=offset,grid=str(grid))
            try:
                prior_path=bundle/'population_prior.json'
                if prior_path.exists():
                    prior=json.loads(prior_path.read_text())
                    if prior['constraint']['source_sha256']!=sha256_file(html) or prior['isochrone_sha256']!=sha256_file(iso) or prior['source_sha256']!=sha256_file(HERE/'empirical_population.py') or prior['quadrature_order']!=order or prior['escape_slope_offset']!=offset:
                        raise ValueError('Frozen empirical population inputs changed')
                else:
                    prior=prepare(cid,iso,reference,html,bundle,offset=offset,order=order)
                update('running_empirical_response_grid',offset=offset,grid=str(grid),population_prior=str(prior_path))
                run_grid(cfg,prior['response_configs'],selection,grid,rates=plan['rates'],times=plan['times_myr'],seeds=plan['seeds'],
                    interval_myr=plan['release_interval_myr'],workers=2,cluster_age_myr=plan['cluster_age_myr'],
                    pixel_factor=plan.get('pixel_factor',1),integrator=plan.get('integrator','cartesian'),batch_size=plan.get('batch_size',32),
                    response_workers=response_workers,reuse_grid=old)
                update('integrating_background_time_population',offset=offset,grid=str(grid))
                result=infer(grid,dest/'posterior',prior_path)
                reports.append(dict(offset=offset,posterior=str(dest/'posterior/posterior_summary.json'),results=result['results']))
                update('offset_completed_conditional_posterior',offset=offset,grid=str(grid),posterior=str(dest/'posterior'))
            except Exception as exc:
                update('empirical_stage_failed',offset=offset,grid=str(grid),error=str(exc));raise
        (root/cid/study_label/'transfer_sensitivity.json').write_text(json.dumps(dict(cluster=cid,
            scenarios=reports,offsets_are_separate_assumptions=True,physical_rate_measurement_validated=False),indent=2),encoding='utf-8')
        update('completed_empirical_posteriors_validation_pending',reports=[r['posterior'] for r in reports])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=HERE/'outputs/cluster_survey_29_v1')
    p.add_argument('--clusters',nargs='+',default=['Pal5','NGC288']);p.add_argument('--offsets',type=float,nargs='+',default=[0.,-.5,-1.])
    p.add_argument('--study-label',default='empirical_mf_v1',help='Use a fresh label for quadrature refinement or revised assumptions')
    p.add_argument('--order',type=int,default=5);p.add_argument('--response-workers',type=int,default=2)
    a=p.parse_args()
    if a.order<3 or a.response_workers<1:p.error('At least three MF nodes and positive worker count required')
    os.environ.setdefault('OPENBLAS_NUM_THREADS','1');run(a.root,a.clusters,a.offsets,a.order,a.response_workers,a.study_label)
