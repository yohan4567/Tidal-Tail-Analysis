"""Prescribed dynamical mass budget, distinct from a self-consistent N-body model."""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class MassEvolution:
    present_mass_msun: float
    duration_myr: float
    stellar_rate_msun_myr: float
    dynamical_to_population_mass: float
    source_description: str

    def __post_init__(self):
        values=[self.present_mass_msun,self.duration_myr,self.stellar_rate_msun_myr,self.dynamical_to_population_mass]
        if not np.isfinite(values).all() or min(values[:2])<=0 or values[2]<0 or values[3]<=0 or not self.source_description.strip():
            raise ValueError('Positive mass/duration/conversion, nonnegative rate and explicit mass-budget assumptions required')

    def mass(self,lookback_myr):
        age=np.asarray(lookback_myr,float)
        if not np.isfinite(age).all() or (age< -1e-7).any() or (age>self.duration_myr+1e-7).any():
            raise ValueError('Mass history evaluated outside [0,duration]')
        age=np.clip(age,0,self.duration_myr)
        return self.present_mass_msun+self.dynamical_to_population_mass*self.stellar_rate_msun_myr*age

    def ratio(self,lookback_myr):
        return self.mass(lookback_myr)/self.present_mass_msun

    def metadata(self):
        return dict(mode='prescribed_linear_escape',present_mass_msun=self.present_mass_msun,
            past_mass_msun=float(self.mass(self.duration_myr)),duration_myr=self.duration_myr,
            stellar_rate_msun_myr=self.stellar_rate_msun_myr,dynamical_to_population_mass=self.dynamical_to_population_mass,
            population_mass_basis='isochrone_current_stellar_mass',amplitude_policy='fixed_reference_rate',
            source_description=self.source_description,plummer_scale_evolution='fixed',
            interpretation='Prescribed M(lookback)=M_now+conversion*stellar_rate*lookback. Release kinematics and post-release forces use the same history.',
            limitations=['Not a relaxation/stellar-evolution prediction of the escape rate.',
                'Conversion between dynamical lost mass and represented present-day living-star mass is assumed.',
                'Plummer scale radius, Galactic potential and progenitor orbit remain prescribed; no dynamical friction or recoil.'])


def from_config(config):
    spec=config.get('mass_evolution')
    if spec is None: return None
    allowed={'mode','dynamical_to_population_mass','source_description','population_mass_basis'}
    if set(spec)-allowed or spec.get('mode')!='linear_escape':
        raise ValueError('Supported mass_evolution is linear_escape with explicitly declared conversion and source')
    if spec.get('population_mass_basis')!='isochrone_current_stellar_mass':
        raise ValueError('Coupled escape requires an explicit current living-stellar-mass basis')
    return MassEvolution(config['mass_msun'],config['duration_myr'],config['reference_mdot_msun_myr'],
        float(spec['dynamical_to_population_mass']),spec['source_description'])
