"""Benchmark matched release states against a stricter scalar reference.

Requires an existing completed simulation; this does not benchmark Gaia inference.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,time
from pathlib import Path
import numpy as np
from dynamics import integrate_released,cylindrical_to_cartesian
from batch_dynamics import integrate_released_batched


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--simulation',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--particles',type=int,default=128);p.add_argument('--layout',choices=['adjacent','spread'],default='adjacent')
    a=p.parse_args()
    if a.particles<2 or a.particles%2:p.error('Positive even particle count required')
    if a.out.exists():p.error('Choose a fresh benchmark output')
    os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).parent/'outputs/.mplcache'))
    import astropy.units as u
    from galpy.orbit import Orbit
    from galpy.potential import MWPotential2014,PlummerPotential
    from galpy.util.conversion import time_in_Gyr
    from mass_evolution import from_config
    from provenance import sha256_file
    meta=json.loads((a.simulation/'run.json').read_text());c=meta['config'];unit=time_in_Gyr(220.,8.)*1000
    with np.load(a.simulation/'release_states.npz') as data:
        count=min(a.particles,len(data['lookback_natural']));half=len(data['lookback_natural'])//2;per=count//2
        start=min(int(.75*half),half-per)
        ids=np.r_[np.arange(start,start+per),np.arange(half+start,half+start+per)] if a.layout=='adjacent' else np.linspace(0,len(data['lookback_natural'])-1,count,dtype=int)
        initial=data['initial_cylindrical'][ids];ages=data['lookback_natural'][ids]
    orbit=Orbit([c['ra_deg'],c['dec_deg'],c['distance_kpc'],c['pmra_masyr'],c['pmdec_masyr'],c['radial_velocity_kms']],radec=True,ro=8.,vo=220.,zo=.0208,solarmotion=[-11.1,12.24,7.25])
    ts=np.linspace(0,-c['duration_myr']/unit,10001)
    orbit.integrate(ts,MWPotential2014,method='dop853_c' if meta['native_acceleration'] else 'dop853',rtol=1e-9,atol=1e-11)
    xyz=np.column_stack([orbit.x(ts,use_physical=False),orbit.y(ts,use_physical=False),orbit.z(ts,use_physical=False)])
    cluster=PlummerPotential(amp=c['mass_msun']*u.Msun,b=meta['plummer_scale_pc']*u.pc,ro=8.,vo=220.)
    evolution=from_config(c)
    amplitude=cluster._amp if evolution is None else lambda t:cluster._amp*evolution.ratio(-np.asarray(t)*unit)
    args=(initial,ages,MWPotential2014,ts,xyz,amplitude,meta['plummer_scale_pc']/8000)
    start=time.perf_counter();cpu=time.process_time();reference=integrate_released(*args,rtol=1e-12,atol=1e-14)
    base=time.perf_counter()-start;base_cpu=time.process_time()-cpu;results={}
    for block in [32,64,128]:
        start=time.perf_counter();cpu=time.process_time();value=integrate_released_batched(*args,batch_size=block,rtol=1e-11,atol=1e-13)
        elapsed=time.perf_counter()-start;cpu_elapsed=time.process_time()-cpu
        delta=cylindrical_to_cartesian(value)-cylindrical_to_cartesian(reference)
        pos=float(np.linalg.norm(delta[:,:3],axis=1).max()*8000);vel=float(np.linalg.norm(delta[:,3:],axis=1).max()*220)
        results[str(block)]=dict(seconds=elapsed,cpu_seconds=cpu_elapsed,speedup=base/elapsed,cpu_speedup=base_cpu/cpu_elapsed,
                                max_position_error_pc=pos,max_velocity_error_kms=vel,passed=pos<.01 and vel<.001)
    a.out.mkdir(parents=True)
    result=dict(tracers=len(ids),layout=a.layout,reference_seconds=base,reference_cpu_seconds=base_cpu,batches=results,
                reference_rtol=1e-12,reference_atol=1e-14,batch_rtol=1e-11,batch_atol=1e-13,
                input_sha256=sha256_file(a.simulation/'release_states.npz'),
                scope='Matched force-law and integration benchmark only; not a physical-rate validation.')
    (a.out/'benchmark.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
