"""Rerun Pal5 with measured tail-MF constraints and verified dynamics reuse.

North/south are separate full-field population-transfer hypotheses, not spatially
matched fits to the original photometric fields, and are never averaged together.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path

from bayesian_grid import run as infer
from provenance import sha256_file
from run_joint_study import run as run_grid
from tail_population import prepare

HERE=Path(__file__).resolve().parent


def run(root,workers=8):
    root=Path(root).resolve();dest=root/'Pal5';source=dest/'production_grid_shared_v1'
    plan=json.loads((source/'plan.json').read_text())
    if not (source/'scores/joint_grid.json').exists():
        raise ValueError('Completed inactive source dynamics required')
    rate,seed=plan['rates'][0],plan['seeds'][0]
    config=source/'configs'/f'rate_{rate:g}_seed_{seed}.json'
    reference=source/'configs'/f'rate_{rate:g}_seed_{seed}_pop_0_T_{max(plan["times_myr"]):g}.json'
    iso=HERE/'outputs/Pal5_science_calibration/isochrone.dat'
    constraint=HERE/'configs/Pal5_tail_mf_Koch2004.json'
    label=dest/'tail_mf_literature_v1';reports=[]
    for region in ('south','north'):
        folder=label/region;bundle=folder/'population_bundle';prior_path=bundle/'population_prior.json'
        grid=folder/'grid'
        if prior_path.exists():
            prior=json.loads(prior_path.read_text())
            if (prior['constraint']['constraint_sha256']!=sha256_file(constraint)
                or prior['constraint']['region']!=region or prior['isochrone_sha256']!=sha256_file(iso)
                or prior['reference_config_sha256']!=sha256_file(reference)
                or prior['source_sha256']!=sha256_file(HERE/'tail_population.py') or prior['quadrature_order']!=5):
                raise ValueError('Frozen literature population changed')
        else:
            prior=prepare(constraint,region,iso,reference,bundle,allow_extrapolation=True)
        def update(stage,**extra):
            (dest/'tail_mf_job.json').write_text(json.dumps(dict(cluster='Pal5',region=region,stage=stage,
                grid=str(grid),workers=workers,updated_utc=datetime.now(timezone.utc).isoformat(),
                total_absolute_rate_supported=False,**extra),indent=2))
        update('running_literature_tail_templates')
        run_grid(config,prior['response_configs'],dest/'selection',grid,rates=plan['rates'],times=plan['times_myr'],
            seeds=plan['seeds'],interval_myr=plan['release_interval_myr'],workers=workers,response_workers=workers,
            cluster_age_myr=plan['cluster_age_myr'],pixel_factor=plan.get('pixel_factor',1),
            integrator=plan['integrator'],batch_size=plan['batch_size'],reuse_grid=source)
        update('integrating_background_time_tail_population')
        infer(grid,folder/'posterior',prior_path)
        report=dict(region=region,posterior=str(folder/'posterior/posterior_summary.json'),
            full_field_population_transfer_assumed=True,spatially_matched_literature_fields=False,
            mass_extrapolation=prior['support_audits'],total_absolute_rate_supported=False)
        (folder/'interpretation.json').write_text(json.dumps(report,indent=2));reports.append(report)
        update('region_completed_conditional_posterior')
    (label/'comparison.json').write_text(json.dumps(dict(reports=reports,
        regions_are_separate_transfer_hypotheses=True,total_absolute_rate_supported=False),indent=2))
    update('completed_literature_tail_posteriors_validation_pending')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=HERE/'outputs/cluster_survey_29_v1')
    p.add_argument('--workers',type=int,default=8);a=p.parse_args()
    if not 1<=a.workers<=8:p.error('Use one to eight workers')
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
    run(a.root,a.workers)
