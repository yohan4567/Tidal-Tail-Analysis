"""Run measured-MF grids directly, omitting legacy tilt grids and posteriors."""
import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path

from empirical_population import prepare, measured_constraint
from production_study import RATES, TIMES, SEEDS
from provenance import sha256_file
from run_joint_study import run as run_grid
from bayesian_grid import run as infer

HERE=Path(__file__).resolve().parent


def run(root, cluster, workers=12, label='empirical_direct_v1'):
    if not 1<=workers<=12:raise ValueError('Use between one and twelve workers')
    root=Path(root).resolve();dest=root/cluster
    if Path(label).name!=label or label in ('.','..'):
        raise ValueError('Study label must be a basename')
    html=HERE/'outputs/empirical_mf_source/parameter.html'
    measured_constraint(html,cluster)  # Fail before starting expensive dynamics.
    rules=json.loads((HERE/'configs/survey_populations_v1.json').read_text())['clusters'].get(cluster,{})
    supplement=HERE/'configs/population_supplements_v1.json'
    if supplement.exists():rules=dict(rules,**json.loads(supplement.read_text())['clusters'].get(cluster,{}))
    if rules.get('hold_reason'):
        raise ValueError(rules['hold_reason'])
    iso=dest/'calibration_inputs/isochrone.dat'
    reference=dest/'science_calibration/bundle/tilt_plus0_zp_minus0.017.json'
    config=dest/'coupled_dynamics.json'
    if not all(p.is_file() for p in (iso,reference,config,dest/'selection/selection.json')):
        raise ValueError('Cluster-specific prepared inputs required')
    reports=[];reuse=dest/'production_grid_shared_v1'
    if not (reuse/'plan.json').exists():reuse=None
    for offset in (0.,-.5,-1.):
        folder=dest/label/('offset_'+f'{offset:g}'.replace('-','minus').replace('.','p'))
        bundle=folder/'population_bundle';prior_path=bundle/'population_prior.json';grid=folder/'grid'
        if prior_path.exists():
            prior=json.loads(prior_path.read_text())
            if (prior['constraint']['source_sha256']!=sha256_file(html) or prior['isochrone_sha256']!=sha256_file(iso)
                or prior['reference_config_sha256']!=sha256_file(reference)
                or prior['source_sha256']!=sha256_file(HERE/'empirical_population.py')
                or prior['quadrature_order']!=5 or prior['escape_slope_offset']!=offset):
                raise ValueError('Frozen population inputs changed')
        else:
            prior=prepare(cluster,iso,reference,html,bundle,offset=offset,order=5)
        def state(stage):
            (dest/'empirical_job.json').write_text(json.dumps(dict(cluster=cluster,stage=stage,offset=offset,
                grid=str(grid),workers=workers,response_workers=workers,
                updated_utc=datetime.now(timezone.utc).isoformat(),physical_rate_measurement_validated=False),indent=2))
        state('running_empirical_response_grid')
        run_grid(config,prior['response_configs'],dest/'selection',grid,rates=RATES,times=TIMES,seeds=SEEDS,
            interval_myr=.25,workers=workers,response_workers=workers,cluster_age_myr=rules['age_gyr']*1000,
            pixel_factor=4,integrator='cartesian_batch',batch_size=128,reuse_grid=reuse)
        state('integrating_background_time_population')
        infer(grid,folder/'posterior',prior_path)
        reports.append(str(folder/'posterior/posterior_summary.json'))
        if offset==0:reuse=grid
    state('completed_empirical_posteriors_validation_pending')
    (dest/label/'transfer_sensitivity.json').write_text(json.dumps(dict(cluster=cluster,reports=reports,
        offsets_are_separate_assumptions=True,physical_rate_measurement_validated=False),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=HERE/'outputs/cluster_survey_29_v1')
    p.add_argument('--cluster',required=True);p.add_argument('--workers',type=int,default=12)
    a=p.parse_args()
    if not 1<=a.workers<=12:p.error('Use between one and twelve workers')
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
    run(a.root,a.cluster,a.workers)
