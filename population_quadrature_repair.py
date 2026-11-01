"""Preserve isochrone phase weights when a rounded mass map has tied nodes."""
import numpy as np


def mass_quadrature(proxy,slope,phase_weights):
    proxy=np.asarray(proxy,float);prior=np.asarray(phase_weights,float)
    if proxy.ndim!=1 or prior.shape!=proxy.shape or len(proxy)<2 or not np.isfinite(proxy).all() or (proxy<=0).any() or (np.diff(proxy)<0).any():
        raise ValueError('A finite positive nondecreasing mass coordinate is required')
    if not np.isfinite(slope) or not np.isfinite(prior).all() or (prior<=0).any():
        raise ValueError('Finite slope and positive original phase weights required')
    unique,inverse=np.unique(proxy,return_inverse=True)
    if len(unique)<2:raise ValueError('Isochrone does not resolve a mass interval')
    exponent=slope+1.;ratio=np.log(unique[1:]/unique[:-1])
    delta=ratio if abs(exponent)<1e-10 else unique[:-1]**exponent*np.expm1(exponent*ratio)/exponent
    group_weights=.5*(np.r_[0.,delta]+np.r_[delta,0.])
    if len(unique)==len(proxy):weights=group_weights
    else:
        totals=np.bincount(inverse,weights=prior,minlength=len(unique))
        weights=group_weights[inverse]*prior/totals[inverse]
    if not np.isfinite(weights).all() or (weights<=0).any():raise ValueError('Invalid tied-mass quadrature')
    if not np.allclose(np.bincount(inverse,weights=weights),group_weights,rtol=2e-14,atol=0):
        raise ValueError('Mass-group measure was not conserved')
    return weights,dict(policy='unique_mass_measure_with_original_phase_allocation_v1',
        tied_nodes=int(len(proxy)-len(unique)),mass_coordinates=int(len(unique)),
        original_nodes=len(proxy),mass_measure=float(group_weights.sum()),
        assumption='Nodes with identical supplied mass proxy share their integrated mass-function measure in proportion to original isochrone phase weights; no photometric node is removed or assigned an invented mass.')
