"""Prespecified evolving-mass pilot grid with matched fixed controls and checkpoints."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
from pathlib import Path

import numpy as np

from provenance import sha256_file
from score_coupled_grid import score
from validate_science import command,check_response

HERE=Path(__file__).resolve().parent


def run(bundle,selection,fixed_simulation,out,*,rates=(4.,9.,25.),particles=2000,workers=2):
    bundle,selection,fixed_simulation,out=[Path(p).resolve() for p in (bundle,selection,fixed_simulation,out)]
    if not np.isfinite(rates).all() or min(rates)<=0 or len(set(rates))!=len(rates) or particles<2 or particles%2:
        raise ValueError('Positive distinct rates and an even tracer count required')
    original=json.loads((fixed_simulation/'run.json').read_text())['config']
    specs=[bundle/f'tilt_{label}_zp_minus0.017.json' for label in ('minus1','plus0','plus1')]
    plan=dict(rates=list(rates),particles=particles,seed=42,source={k:sha256_file(HERE/k) for k in (
        'run_coupled_study.py','score_coupled_grid.py','stream_model.py','dynamics.py','mass_evolution.py','response_model.py','population_model.py','selection_model.py')},
        response_specs={str(p):sha256_file(p) for p in specs},selection_sha256=sha256_file(selection/'selection.json'),
        fixed_simulation_sha256=sha256_file(fixed_simulation/'run.json'),mass_conversion=1.)
    if out.exists():
        if not (out/'plan.json').exists() or json.loads((out/'plan.json').read_text())!=plan:
            raise ValueError('Use a fresh coupled study directory; existing plan differs')
    else:
        out.mkdir(parents=True); (out/'plan.json').write_text(json.dumps(plan,indent=2))
    configs=out/'configs'; configs.mkdir(exist_ok=True)
    simulations={}
    def simulation_job(rate):
        name=f'rate_{rate:g}'; path=out/name/'simulation'; path.parent.mkdir(exist_ok=True)
        c=dict(original,n_particles=particles,seed=42,reference_mdot_msun_myr=rate,
            mass_evolution=dict(mode='linear_escape',dynamical_to_population_mass=1.,population_mass_basis='isochrone_current_stellar_mass',
                source_description='Pilot assumption: all prescribed dynamical mass loss equals represented present-day living-star escape mass; no additional wind/remnant loss; fixed Plummer scale.'))
        cfg=configs/(name+'.json'); cfg.write_text(json.dumps(c,indent=2))
        if not (path/'run.json').exists():
            args=['--config',cfg,'--out',path]
            if path.exists(): args+=['--resume']
            command('stream_model.py',args,out/name/'simulation.log')
        meta=json.loads((path/'run.json').read_text())
        if meta['config']!=dict(c,python_integrator='cartesian'):
            raise ValueError('Simulation mass history mismatch')
        if sha256_file(path/'particles.csv')!=meta['particles_sha256']: raise ValueError('Particle checksum mismatch')
        return rate,path
    with ThreadPoolExecutor(max_workers=workers) as pool:
        tasks=[pool.submit(simulation_job,rate) for rate in rates]
        for f in as_completed(tasks):
            rate,path=f.result(); simulations[rate]=path
            print(f'Completed evolving dynamics: {rate:g} Msun/Myr',flush=True)
    # Common random numbers: changed histories must not change the sampled ages.
    import pandas as pd
    ages=[pd.read_csv(simulations[r]/'particles.csv').lookback_release_myr.to_numpy() for r in rates]
    if not all(np.array_equal(ages[0],a) for a in ages[1:]): raise ValueError('Release sampling changed across histories')
    def response_job(simulation,spec,dest):
        if not (dest/'response.json').exists():
            command('response_model.py',['--selection',selection,'--run',simulation,'--config',spec,'--out',dest],dest.parent/(dest.name+'.log'))
        check_response(dest)
        print('Completed response:',dest.parent.name,dest.name,flush=True)
        return dest/'calibrated_counts.npz'
    coupled=[]; fixed=[]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        tasks={}
        for i,spec in enumerate(specs):
            dest=out/'fixed_controls'/f'population_{i}'
            tasks[pool.submit(response_job,fixed_simulation,spec,dest)]='fixed'
        for rate,path in simulations.items():
            for i,spec in enumerate(specs):
                dest=out/f'rate_{rate:g}'/f'population_{i}'
                tasks[pool.submit(response_job,path,spec,dest)]='coupled'
        for f in as_completed(tasks):
            (fixed if tasks[f]=='fixed' else coupled).append(f.result())
    result=score(coupled,out/'scores',fixed_templates=fixed)
    (out/'REPORT.md').write_text('\n'.join(['# Prescribed evolving-mass pilot grid','',
        'Rates change the release mass and subsequent cluster forces. Only the background is fitted per evolving-rate candidate.',
        '', '![Coupled grid](scores/coupled_grid.png)','',
        'The mass budget uses conversion factor 1 and a fixed Plummer scale. This is a prescribed model, not a self-consistent N-body evolution.',
        'The grid is coarse; no continuous best rate or confidence interval is reported. Population/HB and numerical sampling assumptions remain.',
        '', 'Full candidate scores: scores/scores.csv. Inputs and source hashes: plan.json.'])+'\n',encoding='utf-8')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('bundle','selection','fixed-simulation','out'): p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--rates',type=float,nargs='+',default=[4.,9.,25.]); p.add_argument('--particles',type=int,default=2000)
    p.add_argument('--workers',type=int,default=2); a=p.parse_args()
    r=run(a.bundle,a.selection,a.fixed_simulation,a.out,rates=a.rates,particles=a.particles,workers=a.workers)
    print('Completed coupled pilot:',len(r['candidates']),'scores')
