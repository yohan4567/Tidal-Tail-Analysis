"""Reintegrate matched released states using an independently shortened progenitor orbit."""
import argparse
import json
from pathlib import Path

import numpy as np

from provenance import sha256_file


def validate_window(simulation,duration_myr,max_particles=8,tolerance_pc=.01):
    from galpy.orbit import Orbit
    from galpy.potential import MWPotential2014,PlummerPotential
    from galpy.util.conversion import time_in_Gyr
    import astropy.units as u
    from dynamics import integrate_released
    from mass_evolution import from_config
    path=Path(simulation); meta=json.loads((path/'run.json').read_text()); c=meta['config']
    evolution=from_config(c)
    if evolution is None or not 0<duration_myr<c['duration_myr']:
        raise ValueError('A shorter positive coupled constant-rate window required')
    if max_particles<1 or tolerance_pc<=0: raise ValueError('Positive audit size/tolerance required')
    unit=time_in_Gyr(220.,8.)*1000
    with np.load(path/'release_states.npz',allow_pickle=False) as d:
        age=d['lookback_natural']; initial=d['initial_cylindrical']; final=d['final_cylindrical']
    eligible=np.flatnonzero(age*unit<=duration_myr)
    if not len(eligible): raise ValueError('No matched released states in requested window')
    ids=eligible[np.linspace(0,len(eligible)-1,min(max_particles,len(eligible)),dtype=int)]
    orbit=Orbit([c['ra_deg'],c['dec_deg'],c['distance_kpc'],c['pmra_masyr'],c['pmdec_masyr'],c['radial_velocity_kms']],
        radec=True,ro=8.,vo=220.,zo=.0208,solarmotion=[-11.1,12.24,7.25])
    ts=np.linspace(0,-duration_myr/unit,10001)
    orbit.integrate(ts,MWPotential2014,method='dop853',rtol=1e-10,atol=1e-12)
    xyz=np.column_stack([orbit.x(ts,use_physical=False),orbit.y(ts,use_physical=False),orbit.z(ts,use_physical=False)])
    cluster=PlummerPotential(amp=c['mass_msun']*u.Msun,b=meta['plummer_scale_pc']*u.pc,ro=8.,vo=220.)
    amplitude=lambda t: cluster._amp*evolution.ratio(-np.asarray(t)*unit)
    repeated=integrate_released(initial[ids],age[ids],MWPotential2014,ts,xyz,amplitude,meta['plummer_scale_pc']/8000)
    def positions(states):
        return np.column_stack([states[:,0]*np.cos(states[:,5]),states[:,0]*np.sin(states[:,5]),states[:,3]])
    def velocities(states):
        R,vR,vT,z,vz,phi=states.T
        return np.column_stack([vR*np.cos(phi)-vT*np.sin(phi),vR*np.sin(phi)+vT*np.cos(phi),vz])
    errors=np.linalg.norm(positions(repeated)-positions(final[ids]),axis=1)*8000
    velerrors=np.linalg.norm(velocities(repeated)-velocities(final[ids]),axis=1)*220
    result=dict(duration_myr=duration_myr,audited_particles=len(ids),maximum_position_difference_pc=float(errors.max()),
        maximum_velocity_difference_kms=float(velerrors.max()),position_tolerance_pc=tolerance_pc,
        passed=bool(np.all(errors<tolerance_pc)),release_states_sha256=sha256_file(path/'release_states.npz'),
        scope='Matched release states, independent shortened progenitor interpolation/integration; not a statistical short-run template or population validation.')
    if not result['passed']: raise RuntimeError('Window reintegration audit failed: '+json.dumps(result))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--simulation',type=Path,required=True)
    p.add_argument('--window-myr',type=float,required=True); p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); result=validate_window(a.simulation,a.window_myr)
    if a.out.exists(): raise ValueError('Refusing to overwrite audit')
    a.out.parent.mkdir(parents=True,exist_ok=True); a.out.write_text(json.dumps(result,indent=2))
