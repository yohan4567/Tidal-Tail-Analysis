"""Plan or execute a uniform-resolution Pal 5 follow-up grid.

A production configuration is not a certificate of physical validity.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from run_joint_study import run

RATES=list(map(float,range(2,33)))
TIMES=list(map(float,range(100,2001,100)))
SEEDS=[1009,2027,3037]


def regular_axis(start,stop,step):
    if not np.isfinite([start,stop,step]).all() or start<=0 or stop<start or step<=0:
        raise ValueError('Positive finite axis limits/step required, with stop >= start')
    count=(stop-start)/step
    if not np.isclose(count,round(count),rtol=0,atol=1e-7):
        raise ValueError('Axis range must be an integer multiple of step')
    if count>10000:raise ValueError('Axis exceeds 10001 nodes; revise the explicit plan')
    return (start+step*np.arange(round(count)+1)).tolist()


def assess(scores):
    """Identify grid-edge optima and unresolved template noise without inventing CIs."""
    rows=scores['candidates']; best=min(rows,key=lambda r:r['negative_log_likelihood'])
    rates=sorted({r['rate_msun_myr'] for r in rows})
    times=sorted({r['duration_myr'] for r in rows})
    boundary={key:best[key] in (values[0],values[-1]) for key,values in
              [('rate_msun_myr',rates),('duration_myr',times)]}
    return dict(best_sampled_candidate=best,best_on_boundary=boundary,
        physical_rate_measurement_validated=False,
        status='conditional_grid_requires_validation',
        requirements=['Independent particle-density refinement at the same pixel scale',
                      'Template noise checks including predeclared outer-tail regions',
                      'Profile interval calibration and held-out recovery across nuisance scenarios',
                      'External population mass-function and selection calibration support',
                      'Grid extension or explicit reporting of unresolved boundaries/plateaus'],
        warning='No automatic confidence interval; a boundary optimum is not a measured limit.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('config','selection','out'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--response-configs',type=Path,nargs='+',required=True)
    parser.add_argument('--cluster-age-myr',type=float,required=True)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--confidence95',action='store_true',help='After dynamics/scoring, calibrate matching-grid confidence sets')
    parser.add_argument('--calibration',type=int,default=200)
    parser.add_argument('--validation',type=int,default=300)
    parser.add_argument('--background-scales',type=float,nargs='+',default=[1.])
    parser.add_argument('--execute',action='store_true',help='Run the expensive resumable grid; default prints the plan')
    for name,default in [('rate-min',2.),('rate-max',32.),('rate-step',1.),
                         ('time-min-myr',100.),('time-max-myr',2000.),('time-step-myr',100.)]:
        parser.add_argument('--'+name,type=float,default=default)
    args=parser.parse_args()
    try:
        rates=regular_axis(args.rate_min,args.rate_max,args.rate_step)
        times=regular_axis(args.time_min_myr,args.time_max_myr,args.time_step_myr)
        if any(not np.isclose(t/.25,round(t/.25),atol=1e-7,rtol=0) for t in times):
            raise ValueError('Times must align with the 0.25 Myr release strata')
    except ValueError as exc:parser.error(str(exc))
    if not np.isfinite(args.cluster_age_myr) or args.cluster_age_myr<max(times):
        parser.error('Supplied age bound must cover the requested maximum time')
    particles=2*round(max(times)/.25)
    plan=dict(rates_msun_myr=rates,times_myr=times,seeds=SEEDS,
        release_interval_myr=.25,particles_per_simulation=particles,
        dynamics_runs=len(rates)*len(SEEDS),total_trajectories=particles*len(rates)*len(SEEDS),
        pixel_factor=4,conditional_only=True,grid_nodes=len(rates)*len(times),
        note='T windows reuse each Tmax simulation without renormalizing mass weights. Runtime is machine dependent.')
    if args.confidence95:
        candidates=len(rates)*len(times)*len(args.response_configs)
        experiments=candidates*len(args.background_scales)*(args.calibration+args.validation)
        plan['confidence_calibration_cost']=dict(candidate_nodes=candidates,synthetic_datasets=experiments,
            candidate_likelihood_fits=experiments*candidates,
            caveat='Full per-node calibration can cost substantially more than scoring the observed grid.')
    print(json.dumps(plan,indent=2),flush=True)
    if not args.execute:return
    run(args.config,args.response_configs,args.selection,args.out,rates=rates,times=times,seeds=SEEDS,
        interval_myr=.25,workers=args.workers,cluster_age_myr=args.cluster_age_myr,pixel_factor=4)
    report=json.loads((args.out/'scores/joint_grid.json').read_text())
    (args.out/'production_assessment.json').write_text(json.dumps(assess(report),indent=2),encoding='utf-8')
    if args.confidence95:
        from science_validation import validate_grid
        paths=json.loads((args.out/'templates.json').read_text())['templates']
        validate_grid(paths,args.out/'scores/joint_grid.json',args.out/'validation95',
                      args.calibration,args.validation,background_scales=args.background_scales)


if __name__=='__main__':main()
