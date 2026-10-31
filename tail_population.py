"""Build externally constrained tail populations without assuming a delta.

The output is a conditional population, never a certificate of an absolute rate.
Full-support extrapolation requires an explicit opt-in and is recorded as a blocker.
Existing simulation code and source-frozen grids are not modified.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from empirical_population import measured_population, normal_nodes
from provenance import sha256_file


def read_constraint(path, region):
    data = json.loads(Path(path).read_text(encoding='utf-8'))
    if data.get('convention') != 'dN/dm proportional to m**(-(1+x))':
        raise ValueError('Explicit x convention required; never guess a slope sign')
    if data.get('mass_coordinate') != 'main_sequence_current_mass':
        raise ValueError('Unsupported measured mass coordinate')
    if data.get('measurement_scope') != 'present_day_observed_tail_segment':
        raise ValueError('Expected an observed tail constraint, not a progenitor MF')
    if not data.get('source') or not data.get('cluster'):
        raise ValueError('Cluster and measurement source required')
    lo, hi = data['mass_range_msun']
    r = data['regions'][region]
    if not np.isfinite([lo, hi, r['x'], r['sigma_x']]).all() or not 0 < lo < hi or r['sigma_x'] <= 0:
        raise ValueError('Invalid mass interval or uncertainty')
    return dict(data, region=region, alpha_mean=-(1+r['x']), alpha_sigma=r['sigma_x'],
                constraint_sha256=sha256_file(path))


def support_audit(population, constraint):
    lo, hi = constraint['mass_range_msun']
    p = population.table
    mass = p.current_mass_msun.to_numpy(float)
    weight = p.number_per_msun.to_numpy(float)
    ms = p.phase.to_numpy() == 1
    total = float(np.sum(weight*mass))
    inside = ms & (mass >= lo) & (mass <= hi)
    # A node-weight diagnostic, not a completeness correction or a measured fraction.
    fraction = float(np.sum(weight[inside]*mass[inside])/total)
    minimum, maximum = population.metadata['represented_ms_mass_range_msun']
    extrapolated = minimum < lo-1e-8 or maximum > hi+1e-8
    return dict(measured_mass_range_msun=[lo, hi], represented_ms_mass_range_msun=[minimum, maximum],
                mf_extrapolated=bool(extrapolated), model_mass_fraction_at_measured_ms_nodes=fraction,
                evolved_phase_transfer_assumed=bool(np.any(~ms)),
                total_absolute_rate_supported=False,
                rate_scope=population.metadata['rate_scope'])


def prepare(constraint_path, region, isochrone, reference_config, out, *, order=5,
            allow_extrapolation=False):
    constraint_path, isochrone, reference_config, out = map(Path, (constraint_path, isochrone, reference_config, out))
    if out.exists():
        raise ValueError('Use a fresh output directory; existing grids are immutable')
    constraint = read_constraint(constraint_path, region)
    imeta = json.loads(isochrone.with_suffix('.json').read_text(encoding='utf-8'))
    if imeta['sha256'] != sha256_file(isochrone):
        raise ValueError('Isochrone checksum changed')
    spec = json.loads(reference_config.read_text(encoding='utf-8'))
    base = reference_config.resolve().parent
    for key in ('quality', 'cache_dir', 'rate_history'):
        if spec.get(key):
            spec[key] = str((base/spec[key]).resolve())
    for section in ('errors', 'completeness'):
        for key in ('path', 'metadata'):
            if isinstance(spec[section].get(key), str):
                spec[section][key] = str((base/spec[section][key]).resolve())
    slopes, weights = normal_nodes(constraint['alpha_mean'], constraint['alpha_sigma'], order)
    populations, audits = [], []
    for slope in slopes:
        pop = measured_population(isochrone, slope, age=imeta['age_years'], mh=imeta['metallicity_mh'],
                                  source=f"{constraint['citation']}; {region}; directly measured tail MF")
        audit = support_audit(pop, constraint)
        if audit['mf_extrapolated'] and not allow_extrapolation:
            raise ValueError('Measured tail MF does not cover represented masses; explicit --allow-extrapolation only creates a conditional model')
        # Replace inherited progenitor-transfer wording: this route does not use delta.
        assumptions = [a for a in pop.metadata['assumptions']
                       if not a.startswith(('Single global-MF', 'Escape MF equality/offset'))]
        assumptions += ['Measured tail-segment slope is extended over the full represented mass range.',
                        'Observed segment MF is assumed representative of this whole model; transfer unvalidated.',
                        'Published formal slope error is approximated by a Gaussian; systematic uncertainty is not included.']
        pop.metadata.update(mass_function_status='Measured present-day tail-segment slope; conditional mass/spatial/time extrapolation',
                            tail_mf_constraint=constraint, support_audit=audit, assumptions=assumptions,
                            tail_population_source_sha256=sha256_file(Path(__file__)))
        populations.append(pop)
        audits.append(audit)
    out.mkdir(parents=True)
    configs, prior_weights = [], {}
    for i, (pop, slope, weight) in enumerate(zip(populations, slopes, weights)):
        pop.save(out/f'population_{i}')
        scenario = dict(measured_tail_mf_slope=float(slope), tail_region=region,
                        tail_constraint_sha256=constraint['constraint_sha256'],
                        parallax_zero_point_mas=spec['parallax_zero_point_mas'])
        response = dict(spec, population=str((out/f'population_{i}').resolve()), scenario=scenario,
                        source_description=f"{constraint['citation']}; {region}; measured tail MF with explicit conditional extrapolation")
        response.pop('lookback_window_myr', None)
        path = out/f'response_{i}.json'
        path.write_text(json.dumps(response, indent=2), encoding='utf-8')
        configs.append(str(path.resolve()))
        prior_weights[json.dumps(scenario, sort_keys=True)] = float(weight)
    report = dict(cluster=constraint['cluster'], constraint=constraint, quadrature_order=order,
                  scenario_weights=prior_weights, response_configs=configs,
                  support_audits=audits, total_absolute_rate_supported=False,
                  measurement_reuse='Use the published tail MF once as external information; do not also multiply the same photometry likelihood.',
                  population_model_scope='Represented living stars; excludes remnants, winds and stars outside the isochrone.',
                  blockers=['Unmeasured MF outside fitted mass interval', 'Tail-segment to full escape-history transfer',
                            'Remnant and low-mass closure', 'Independent selection/background calibration',
                            'Dynamics, numerical convergence and injection recovery'],
                  source_sha256=sha256_file(Path(__file__)), isochrone_sha256=sha256_file(isochrone),
                  reference_config_sha256=sha256_file(reference_config))
    (out/'population_prior.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('constraint', 'isochrone', 'reference-config', 'out'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--region', required=True)
    p.add_argument('--order', type=int, default=5)
    p.add_argument('--allow-extrapolation', action='store_true')
    a = p.parse_args()
    result = prepare(a.constraint, a.region, a.isochrone, a.reference_config, a.out,
                     order=a.order, allow_extrapolation=a.allow_extrapolation)
    print(json.dumps(dict(out=str(a.out), total_absolute_rate_supported=False, blockers=result['blockers']), indent=2))
