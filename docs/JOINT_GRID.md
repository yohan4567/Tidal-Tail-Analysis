# Conditional joint time and mass loss grid

This implements the first constant-rate stage of [the version 3 method](EXPERIMENTAL_METHOD_ZH.md).
It does not certify a physical mass-loss measurement. Existing fixed-dynamics outputs remain separate.

## Run

Use PowerShell from `gc_stream_model`. The dynamics JSON must explicitly contain
`mass_evolution` with `mode: linear_escape`, the population mass basis, conversion factor and source.
Use an independently supported cluster age and assess its uncertainty and environmental history.
The example age is the existing conditional Pal5 population assumption, not a new age measurement.

```powershell
.\.venv\Scripts\python.exe run_joint_study.py --config outputs/evolving_mass_smoke_config.json --response-configs outputs/Pal5_science_calibration/bundle/tilt_minus1_zp_minus0.017.json outputs/Pal5_science_calibration/bundle/tilt_plus0_zp_minus0.017.json outputs/Pal5_science_calibration/bundle/tilt_plus1_zp_minus0.017.json --selection outputs/Pal5_real_inputs_baseline/selection --out outputs/Pal5_joint_constant_v1 --rates 4 9 25 --times-myr 500 1000 2000 --interval-myr 1 --seeds 42 84 --cluster-age-myr 11000 --workers 2
```

This is a computational experiment with conditional mass conversion and radius assumptions.
Review `applicability.json`; large past/current mass ratios are warnings, not automatic physical exclusions.
Choose a new output folder for a changed plan. An unchanged plan can resume completed steps;
partial dynamics use verified checkpoints. An incomplete response directory must be inspected and
recovered separately; it is not silently overwritten.

## Sampling and reuse

There is one Tmax integration per rate and independent dynamics seed, not per T/rate cell.
`release_sampling: stratified` draws one time within each equal-width bin per arm.
The driver sets N = 2 Tmax / interval. Every requested time must end on a bin boundary.
Each tracer has original mass weight rate * Tmax / N; shorter windows set older weights to zero.
No subset renormalization is performed. This is a modeled escape-population window, not a formation event.

The sampler uses galpy's inverse-CDF hook separately per arm; this is a stratified Monte Carlo design,
not a claim about bursts of physical escape. `uniform` remains supported for existing runs.

The first dynamics realization undergoes a matched-state short-window reintegration audit,
using an independently computed shorter progenitor interpolation. `window_audit.json` records the
position and velocity differences. This audits a small set of trajectories; it does not replace
independently sampled short-run template/refinement validation. Other families or histories require
their own audits. Source hashes and seeds remain in the plan.

## Observation response and age preflight

`response_model.py` accepts `lookback_window_myr` in its JSON. Shortening a window currently requires
a prescribed linear coupled history. Arbitrary histories cannot accidentally pass this guard.
It writes `age_contributions.csv`: selected fitted-region counts by release-age bin and cumulative
older contribution. `tracer_response.npz` preserves finer tracer-age contributions if different age
bins are needed. Weights include the declared population, survey/quality and science selection response.
These are predicted counts, not observed age assignments.

The driver now computes the full observation response only once per dynamics/population pair.
`response_windows.py` sums its saved per-tracer, per-calibration contributions to make shorter T
templates. It checks that these contributions reconstruct the full fitted-region template before
reusing them. Thus extra T points cost count aggregation and scoring, rather than repeated observation
draws. Derived maps are fitted-region maps; outside-fit signal values are zero. Their per-observation
MC error was not stored and is explicitly unavailable (NaN), rather than copied from the full window.

Use a Tmax response on representative rates/populations as a cheap preflight before a large grid.
Inspect cumulative older contributions and changes in whole-map predictions; a T plateau does not
automatically imply a lower bound. The driver currently generates all requested windows rather than
automatically choosing or truncating the T grid.

## Scores and template noise

`joint_grid.py --templates ... --out ...` profiles background only at unit signal amplitude.
It writes the sampled T/rate grid, rate/time profiles, background-only score and an image.
Population scenarios are profiled as a discrete declared set. The data, background, mass basis,
model version and dynamics family must agree. There is no interpolated optimum or calibrated CI.

The background fit uses `joint_fit.py`. A mixture of shared Poisson calibration states can have
multiple maxima even though each component is log-concave. The solver bounds the mixture using
component maxima and concave tangents, subdividing until the log-likelihood bound is within 1e-8.
It raises an error rather than claiming convergence if its iteration/precision limit is reached.
Single bounded scalar minimization is not sufficient for arbitrary mixtures.

Multiple independent dynamics seeds are averaged with aligned response/calibration states. Their
mean-template variation is diagnosed separately from calibration-state uncertainty. At least two
distinct dynamics seeds are needed for a diagnostic; two seeds alone do not certify convergence.
The top-90%-signal pixel target is SE(template mean) < 0.1 sqrt(expected signal + background).
Pass a predeclared `--regions regions.npz` with named boolean masks in fitted-pixel order to also
check distant-tail regional sums while retaining spatial covariance. Masks must not be chosen from
the observed residuals. No kernel smoothing is applied automatically.

An independent density/refinement comparison and stability of likelihood/profile results remain
necessary. Repeated observation draws of the same trajectory are never called independent dynamics
particles. Results always remain conditional; no automated promotion to reporting level 1 exists.

## Injection recovery

`joint_recovery.py --templates ... --manifest injections.json --out ... --trials 100` accepts a manifest:

```json
{
  "source_description": "Predeclared frozen-template experiment",
  "scenarios": [
    {"name": "interior", "template": "path/relative/to/manifest/template.npz", "background_scale": 1, "tags": ["interior"]},
    {"name": "background_only", "null": true, "tags": ["null"]}
  ]
}
```

Also supply boundary, low_signal, time_sensitive, time_transition, time_saturated and off_grid
scenarios; missing tags are reported. Tags are declarations and must be justified by preflight.
Independent truth templates are preferred. A frozen-template test does not rerun selection or
calibration construction and is explicitly labeled accordingly.

Trials are split into calibration and held-out halves. For on-grid truth the code calibrates a
scenario-specific empirical 95% q threshold on one half and evaluates point inclusion on the other.
For a null case it similarly checks a held-out false-positive fraction. Off-grid cases report bias,
not unsupported continuous-region coverage. Thresholds cannot be reused as universal observed-data
confidence cutoffs. At least 40 trials are enforced; more are needed for precise coverage assessment.
The reported binomial standard errors are approximate, particularly near zero or one; Wilson
intervals are also supplied. Recovery uses the same aligned, seed-averaged fitted templates as the
observed-data analysis. It no longer silently uses only the first dynamics seed.

## Executed validation

`validate_joint_science.py` runs the prespecified Pal5 experiment and saves sources, inputs and
tolerances. `finalize_joint_validation.py` can rescore frozen completed response stages with the
globally bounded profiler, preserving the raw products and recording the analysis source hashes.
The first raw experiment's original scores are retained for provenance and superseded by the
corrected final analysis, not silently overwritten.

The executed validation uses 3 rates, 4 windows, 3 population assumptions and 3 dynamics seeds;
it includes doubled particle density, independently generated short-window dynamics, an off-grid
physical rate and frozen-template injection recovery. Passing software checks does not imply that
the finite-particle spatial templates pass the stricter convergence target. See the actual final
report under `outputs/Pal5_joint_validation_final/REPORT_ZH.md` before interpreting a rate.

## Remaining scientific validation

This code supplies tools, not completed research evidence: choose predeclared regions/tolerances,
run adequate independent/refined dynamics, review age preflight and applicability, populate the
injection matrix, test population assumptions, and validate against independent dynamics where possible.
No self-consistent radius/relaxation evolution, automatically calibrated population mass function,
or time-varying escape-rate inference is claimed.
