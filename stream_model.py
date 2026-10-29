"""Chen et al. stream dynamics; tracer mass is not a Gaia star count."""
import argparse
import csv
import hashlib
import json
import os
import platform
import time
from pathlib import Path

import numpy as np
from mass_evolution import from_config
from provenance import sha256_file

HERE = Path(__file__).resolve().parent


def validate(c):
    for key in ('distance_kpc', 'mass_msun', 'half_mass_radius_pc', 'duration_myr',
                'reference_mdot_msun_myr'):
        if not np.isfinite(c[key]) or c[key] <= 0:
            raise ValueError(f'{key} must be positive and finite')
    for key in ('ra_deg', 'dec_deg', 'pmra_masyr', 'pmdec_masyr', 'radial_velocity_kms'):
        if not np.isfinite(c[key]):
            raise ValueError(f'{key} must be finite')
    if not 0 <= c['ra_deg'] < 360 or not -90 <= c['dec_deg'] <= 90:
        raise ValueError('Invalid RA/Dec')
    for key in ('n_particles',):
        if isinstance(c[key], bool) or int(c[key]) != c[key] or c[key] < 1:
            raise ValueError(f'{key} must be a positive integer')
    if c['n_particles'] < 2 or c['n_particles'] % 2:
        raise ValueError('n_particles must be even, at least 2 (equal arms)')
    if c['potential'] != 'MWPotential2014':
        raise ValueError('This version supports MWPotential2014 only')
    from_config(c)
    if c.get('release_sampling','uniform') not in ('uniform','stratified'):
        raise ValueError('release_sampling must be uniform or stratified')


def mass_weights(n, duration, mdot):
    """Uniform release-time sampling; total mass covers BOTH arms."""
    if n < 1 or duration <= 0 or mdot < 0:
        raise ValueError('Invalid mass normalization')
    return np.full(n, duration * mdot / n)


def write_columns(path, columns):
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(columns)
        w.writerows(zip(*columns.values()))


def simulate(c, out, resume=False):
    started=time.perf_counter()
    os.environ.setdefault('MPLCONFIGDIR',str(HERE/'outputs/.mplcache'))
    import astropy
    import astropy.units as u
    import astropy.coordinates as ac
    import galpy
    from galpy.orbit import Orbit
    from galpy.potential import MWPotential2014, PlummerPotential,TimeDependentAmplitudeWrapperPotential
    from galpy.df import chen24spraydf
    from galpy.util.conversion import time_in_Gyr
    from galpy.util._load_extension_libs import load_libgalpy
    from unittest.mock import patch
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    validate(c)
    evolution=from_config(c)
    out = Path(out)
    signature=dict(config=c,source_sha256={name:sha256_file(HERE/name) for name in ('stream_model.py','dynamics.py','mass_evolution.py')},
        versions=dict(numpy=np.__version__,galpy=galpy.__version__,astropy=astropy.__version__))
    if resume:
        if not (out/'simulation_input.json').exists() or json.loads((out/'simulation_input.json').read_text())!=signature:
            raise ValueError('Resume requires identical config, source and versions in simulation_input.json')
    else:
        out.mkdir(parents=True, exist_ok=False)
        (out/'simulation_input.json').write_text(json.dumps(signature,indent=2),encoding='utf-8')
    # Use the native physical normalization of MWPotential2014 (8 kpc, 220 km/s).
    ro, vo = 8.0, 220.0
    solar = [-11.1,12.24,7.25]  # galpy convention: first component points outward
    sky = ac.SkyCoord(ra=c['ra_deg']*u.deg, dec=c['dec_deg']*u.deg,
                      distance=c['distance_kpc']*u.kpc,
                      pm_ra_cosdec=c['pmra_masyr']*u.mas/u.yr,
                      pm_dec=c['pmdec_masyr']*u.mas/u.yr,
                      radial_velocity=c['radial_velocity_kms']*u.km/u.s)
    prog = Orbit([c['ra_deg'],c['dec_deg'],c['distance_kpc'],c['pmra_masyr'],
                  c['pmdec_masyr'],c['radial_velocity_kms']],radec=True,
                 ro=ro,vo=vo,zo=0.0208,solarmotion=solar)
    mw = MWPotential2014
    a_pc = c['half_mass_radius_pc'] * np.sqrt(2**(2/3)-1)
    cluster = PlummerPotential(amp=c['mass_msun']*u.Msun,b=a_pc*u.pc,ro=ro,vo=vo)
    time_unit_myr=time_in_Gyr(vo,ro)*1000
    mass_argument=c['mass_msun']*u.Msun
    cluster_for_sampling=cluster
    amplitude=cluster._amp
    if evolution is not None:
        amplitude=lambda t: float(cluster._amp*evolution.ratio(-np.asarray(t)*time_unit_myr)) if np.ndim(t)==0 else cluster._amp*evolution.ratio(-np.asarray(t)*time_unit_myr)
        mass_argument=amplitude  # galpy callable with natural-unit input/output.
        cluster_for_sampling=TimeDependentAmplitudeWrapperPotential(pot=cluster,A=lambda t: float(evolution.ratio(-t*time_unit_myr)),ro=ro,vo=vo)
    np.random.seed(c['seed'])  # galpy's sampler uses numpy's global RNG
    print(f"Generating {c['name']}: {c['n_particles']} tracers, Chen DF, {c['potential']}", flush=True)
    native=bool(load_libgalpy()[1])
    method='dop853_c' if native else 'dop853'
    original_integrate=Orbit.integrate
    original_call=Orbit._call_internal
    def call_with_scalar_time(self,*args,**kwargs):
        # galpy 1.12's moving potential passes a dimensionless 0-D array.
        # Preserve quantities and vector times; only unwrap scalar ndarrays.
        if args and type(args[0]) is np.ndarray and args[0].ndim == 0:
            args=(args[0].item(),*args[1:])
        if type(kwargs.get('t')) is np.ndarray and kwargs['t'].ndim == 0:
            kwargs['t']=kwargs['t'].item()
        return original_call(self,*args,**kwargs)
    def integrate_with_method(self,t,pot,*args,**kwargs):
        # The upstream sampler does not expose its internal integration method.
        # This scoped override changes only numerical integration, not the Chen DF.
        kwargs['method']=method
        kwargs.setdefault('rtol',1e-9); kwargs.setdefault('atol',1e-11)
        kwargs.setdefault('progressbar',False)
        return original_integrate(self,t,pot,*args,**kwargs)
    print(f'Integrator: {method} (native acceleration: {native})',flush=True)
    with patch.object(Orbit,'integrate',integrate_with_method), patch.object(Orbit,'_call_internal',call_with_scalar_time):
        df = chen24spraydf(mass_argument,progenitor=prog,pot=mw,progpot=cluster_for_sampling,
                          tdisrupt=c['duration_myr']*u.Myr,tail='both',ro=ro,vo=vo)
        if evolution is not None and not np.isclose(df._progenitor_mass_fn(-df._tdisrupt),amplitude(-df._tdisrupt)):
            raise RuntimeError('Installed galpy does not apply the prescribed release mass history')
        if c.get('release_sampling')=='stratified':
            # galpy calls this once per arm, passing independent uniforms.
            # Each arm samples one time in every equal-width stratum.
            if not hasattr(df,'_stripping_inv_cdf'):
                raise RuntimeError('Installed galpy lacks release-time inverse CDF support')
            df._stripping_inv_cdf=lambda u: -(np.arange(len(u))+u)/len(u)*df._tdisrupt
        if native and evolution is None:
            stream, dt = df.sample(c['n_particles'],returndt=True)
        else:
            # scipy's bounded solver never evaluates beyond the progenitor's
            # interpolation interval (galpy's Python dop853 can overshoot it).
            from scipy.integrate import solve_ivp
            from galpy.potential import evaluateRforces, evaluatezforces, evaluatephitorques
            released, dt = df.sample(c['n_particles'],returndt=True,integrate=False)
            states=[]
            def rhs(t,y):
                R,vR,vT,z,vz,phi=y
                return [vR,vT*vT/R+evaluateRforces(df._pot,R,z,phi=phi,t=t,use_physical=False),
                        -vR*vT/R+evaluatephitorques(df._pot,R,z,phi=phi,t=t,use_physical=False)/R,
                        vz,evaluatezforces(df._pot,R,z,phi=phi,t=t,use_physical=False),vT/R]
            if c.get('python_integrator','cartesian')=='cartesian':
                from dynamics import integrate_released
                ts=df._progenitor_times
                xyz=np.column_stack([df._progenitor.x(ts,use_physical=False),df._progenitor.y(ts,use_physical=False),df._progenitor.z(ts,use_physical=False)])
                states=integrate_released(released.vxvv,dt,mw,ts,xyz,amplitude,a_pc/(1000*ro),
                    checkpoint_path=out/'integration_checkpoint.npz',model_identity=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest())
                method='scipy_DOP853_cartesian_cached_trajectory'
            elif c.get('python_integrator')=='galpy_reference':
                for i,(state,age) in enumerate(zip(released.vxvv,dt)):
                    result=solve_ivp(rhs,(-float(age),0.),state,method='DOP853',rtol=1e-9,atol=1e-11)
                    if not result.success:
                        raise RuntimeError(result.message)
                    states.append(result.y[:,-1])
                    if (i+1)%20==0:
                        print(f'Integrated {i+1}/{len(dt)} tracers',flush=True)
                method='scipy_DOP853 (particles); galpy_dop853 (progenitor)'
            else:
                raise ValueError('Unknown Python integrator')
            stream=Orbit(np.asarray(states),ro=ro,vo=vo,zo=0.0208,solarmotion=solar)
            np.savez_compressed(out/'release_states.npz',initial_cylindrical= released.vxvv,
                lookback_natural=np.asarray(dt),final_cylindrical=np.asarray(states))
    # With return_orbit=True galpy returns dt in its internal time units.
    lookback = np.asarray(dt)*time_in_Gyr(vo,ro)*1000
    weights = mass_weights(c['n_particles'],c['duration_myr'],c['reference_mdot_msun_myr'])
    observed = stream.SkyCoord()
    offset = observed.transform_to(ac.SkyOffsetFrame(origin=ac.SkyCoord(ra=sky.ra, dec=sky.dec)))
    columns = dict(ra_deg=observed.ra.deg, dec_deg=observed.dec.deg,
                   distance_kpc=observed.distance.to_value(u.kpc),
                   pmra_masyr=observed.pm_ra_cosdec.to_value(u.mas/u.yr),
                   pmdec_masyr=observed.pm_dec.to_value(u.mas/u.yr),
                   radial_velocity_kms=observed.radial_velocity.to_value(u.km/u.s),
                   lookback_release_myr=lookback, mass_weight_msun=weights,
                   arm=np.array(['leading']*(len(weights)//2)+['trailing']*(len(weights)//2)),
                   offset_lon_deg=offset.lon.wrap_at(180*u.deg).deg, offset_lat_deg=offset.lat.deg)
    if evolution is not None:
        columns['progenitor_mass_at_release_msun']=evolution.mass(lookback)
        grid=np.linspace(0,c['duration_myr'],1001)
        write_columns(out/'mass_history.csv',dict(lookback_myr=grid,mass_msun=evolution.mass(grid)))
    if not all(np.all(np.isfinite(v)) for k,v in columns.items() if k!='arm'):
        raise RuntimeError('Nonfinite simulation output')
    write_columns(out/'particles.csv', columns)
    times=np.linspace(0,-c['duration_myr']/1000,1001)*u.Gyr
    orbit_method='dop853_c' if native else 'dop853'
    prog.integrate(times,mw,method=orbit_method,rtol=1e-9,atol=1e-11)
    op=np.array([prog.x(times),prog.y(times),prog.z(times)])
    write_columns(out/'orbit.csv', dict(time_myr=times.to_value(u.Myr),
                                      x_kpc=op[0], y_kpc=op[1], z_kpc=op[2]))
    back=prog(times[-1])
    back.integrate(np.linspace(0,c['duration_myr']/1000,1001)*u.Gyr,mw,method=orbit_method,rtol=1e-9,atol=1e-11)
    final=c['duration_myr']/1000*u.Gyr
    position_error=float(np.linalg.norm([back.x(final)-prog.x(0*u.Gyr),back.y(final)-prog.y(0*u.Gyr),back.z(final)-prog.z(0*u.Gyr)])*1000)
    velocity_error=float(np.linalg.norm([back.vx(final)-prog.vx(0*u.Gyr),back.vy(final)-prog.vy(0*u.Gyr),back.vz(final)-prog.vz(0*u.Gyr)]))
    # Verify normalization against the actual sampled release interval.
    represented_time = c['duration_myr']
    expected_mass = c['reference_mdot_msun_myr']*represented_time
    if not np.isclose(weights.sum(), expected_mass):
        raise RuntimeError('Tracer mass normalization failed')
    meta = dict(config=c, integrator=method,native_acceleration=native,plummer_scale_pc=a_pc, n_tracers=len(weights),
                model_version='prescribed_mass_history_v1' if evolution is not None else 'fixed_dynamics_amplitude_v1',
                dynamics_mass_history=None if evolution is None else evolution.metadata(),
                release_time_sampling='stratified_lookback_both_arms' if c.get('release_sampling')=='stratified' else 'uniform_lookback_both_arms',
                represented_duration_myr=represented_time, represented_mass_msun=float(weights.sum()),
                mass_fraction_vs_fixed_progenitor=float(weights.sum()/c['mass_msun']),
                present_position_error_pc=position_error, present_velocity_error_kms=velocity_error,
                solar_frame=dict(ro_kpc=ro,vo_kms=vo,z_sun_pc=20.8,galpy_solar_peculiar_kms=solar),
                versions=dict(python=platform.python_version(),numpy=np.__version__,astropy=astropy.__version__,galpy=galpy.__version__),
                interpretation='Forward-model tracer mass normalization, NOT an observed mass-loss measurement.',
                limitations=['Fixed progenitor mass and Plummer radius; no self-consistent mass evolution.',
                             'No Gaia photometry, selection completeness or measurement-error model.',
                             'No tracer-tracer gravity or internal collisional relaxation.',
                             'Large represented mass fractions require evolving progenitor models.'])
    meta['config_sha256'] = hashlib.sha256(json.dumps(c,sort_keys=True).encode()).hexdigest()
    meta['source_sha256']=signature['source_sha256']
    meta['particles_sha256']=sha256_file(out/'particles.csv')
    if evolution is not None:
        meta['limitations'][0]='Prescribed evolving progenitor mass; fixed Plummer radius; no internal relaxation or stellar-evolution feedback.'
        meta['limitations'][-1]='Each physical reference rate requires its own dynamics; arbitrary amplitude rescaling is not a coupled rate inference.'
    (out/'run.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    fig, axes = plt.subplots(1,3,figsize=(15,4.5),constrained_layout=True)
    good = weights>0
    sc=axes[0].scatter(columns['offset_lon_deg'][good],columns['offset_lat_deg'][good],
                       c=columns['lookback_release_myr'][good],s=3,cmap='viridis')
    axes[0].scatter(0,0,marker='*',color='red',s=80)
    axes[0].set(xlabel='Sky offset longitude (deg)',ylabel='Sky offset latitude (deg)',title='Chen stream (no Gaia selection)')
    fig.colorbar(sc,ax=axes[0],label='Release lookback time (Myr)')
    axes[1].scatter(columns['pmra_masyr'][good],columns['pmdec_masyr'][good],s=3)
    axes[1].scatter(c['pmra_masyr'],c['pmdec_masyr'],marker='*',color='red',s=80)
    axes[1].set(xlabel='PM RA*cos(dec) (mas/yr)',ylabel='PM Dec (mas/yr)',title='Intrinsic proper motions')
    axes[2].plot(op[0],op[1]); axes[2].scatter(op[0,0],op[1,0],color='red')
    axes[2].set(xlabel='Galactocentric x (kpc)',ylabel='Galactocentric y (kpc)',title='Progenitor orbit',aspect='equal')
    axes[2].set_aspect('equal',adjustable='datalim')
    fig.suptitle(f"{c['name']} | assumed rate {c['reference_mdot_msun_myr']:g} Msun/Myr | {'prescribed evolving mass' if evolution else 'fixed mass'}")
    fig.savefig(out/'overview.png',dpi=160); plt.close(fig)
    meta['elapsed_seconds']=time.perf_counter()-started
    (out/'run.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    print(f'Output: {out.resolve()}\nTracers: {len(weights)}; represented mass: {weights.sum():.3f} Msun',flush=True)
    return meta


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,default=HERE/'configs/Pal5.json')
    p.add_argument('--out',type=Path,default=HERE/'outputs/Pal5')
    p.add_argument('--duration-myr',type=float)
    p.add_argument('--particles',type=int)
    p.add_argument('--mdot',type=float,help='Assumed total rate over BOTH arms (Msun/Myr)')
    p.add_argument('--python-integrator',choices=('cartesian','galpy_reference'),default='cartesian')
    p.add_argument('--resume',action='store_true',help='Resume a checkpoint with identical config/source/versions; stop the old process first')
    args=p.parse_args()
    c=json.loads(args.config.read_text(encoding='utf-8'))
    if args.duration_myr is not None: c['duration_myr']=args.duration_myr
    if args.particles is not None: c['n_particles']=args.particles
    if args.mdot is not None: c['reference_mdot_msun_myr']=args.mdot
    c['python_integrator']=args.python_integrator
    simulate(c,args.out,resume=args.resume)


if __name__=='__main__':
    main()
