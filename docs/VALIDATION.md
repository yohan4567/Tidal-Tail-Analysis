# Conditional scientific validation

Run after preparing the real-input calibration bundle and baseline analysis:

```powershell
.\.venv\Scripts\python.exe validate_science.py --bundle outputs/Pal5_science_calibration/bundle --baseline outputs/Pal5_real_inputs_baseline --simulation outputs/Pal5_hybrid_v2/simulation --out outputs/Pal5_science_validation --workers 2 --particles 4000 --gof-trials 200
```

The suite saves its plan (input/configuration/source hashes and thresholds)
before executing. Two independent local Python processes can run at once;
each writes to its own directory. It reuses the existing baseline and checks
completed response hashes on resume. Rerun the identical command after an
interruption; changed inputs, code or a different plan require a new output
directory. Do not manually replace a fitted template or result.
Stop the previous validation processes before resuming. An empty interrupted
simulation placeholder can be recovered automatically; a directory containing
partial scientific data is preserved and requires a fresh output directory.

## What it checks

1. **Nine assumed-model scenarios.** The escape-population weight tilts -1/0/+1
   are crossed with scalar parallax offsets -0.030/-0.017/0 mas. Every scenario
   reruns the observation response and fits both the signal amplitude and
   background normalization on the identical observed count map. The observed
   sample, background shape, dynamics and extinction are held fixed. Shared
   calibration states are averaged at the entire-data likelihood level.

   The rate range across these scenarios is not a statistical confidence
   interval. Similar likelihood values do not establish a mass function: its
   normalization can be almost degenerate with the escaped-mass amplitude.
   Scalar offsets are sensitivity settings, not the full per-source Gaia
   zero-point correction.

2. **Observation Monte Carlo resolution.** The baseline's observation draws
   are doubled (32 to 64 for the supplied bundle); its 8 shared quality-map
   draws remain fixed. The fitted rate, total selected signal and normalized
   spatial template are compared. This checks numerical integration of the
   observation response, not all calibration uncertainties.

3. **Tracer resolution.** A new 4,000-tracer simulation (seed 84) is compared
   with the saved 2,000-tracer baseline (seed 42), followed by the same
   observation model and fit. This includes stochastic tracer sampling, not
   only integration accuracy. With just these two realizations, agreement is
   a diagnostic rather than an established asymptotic convergence sequence.

The suite records a prespecified 5% tolerance on rate and integrated selected
signal. It also reports normalized-map total variation, half the sum of the
absolute difference between normalized pixel templates. A small integrated
count change can coexist with a large change in the spatial distribution.
Consequently the boolean tolerance result must not be read as validating all
spatial bins or all aspects of the model.

## Conditional goodness of fit

```powershell
.\.venv\Scripts\python.exe fit_diagnostics.py --input outputs/Pal5_real_inputs_baseline/response/calibrated_counts.npz --out outputs/Pal5_conditional_gof --trials 200
```

The statistic is twice the difference between the saturated-count log likelihood
and the fitted calibration-mixture log likelihood. For a single calibration
state it reduces to Poisson deviance. Factorial constants cancel within each
dataset. Zero-count cells are included.

At the fitted amplitudes, each realization draws **one shared calibration
state for the whole map**, then independent Poisson counts conditional on that
state. Both amplitudes are refitted with the same mixture likelihood. The
fraction with deviance at least as large as observed is reported with a
plus-one correction `(k+1)/(B+1)`, its finite simulation resolution and an exact
binomial interval for the underlying exceedance probability.

**This is a conditional discrepancy check, not a full pipeline p-value.**
The background shape was learned from this field and is held fixed in these
mocks. The selection, empirical CMD ridge, masks, binning and search decisions
are also fixed. A full calibration must regenerate data and repeat those
estimation/selection/search stages. No asymptotic chi-square tail probability
is computed for the sparse Pearson statistic. Ordinary R-squared is saved
only as a descriptive residual statistic.

## Outputs and interpretation

`REPORT.md`, `validation.json`, `sensitivity.csv` and `sensitivity.png` summarize
the run. `models/` holds the individual responses and conditional fits;
`numerics/` holds the observation/tracer comparisons. `goodness_of_fit/` holds
the bootstrap statistics, residual counts, histogram and limitations.

Neither a non-small conditional p-value nor numerical agreement establishes
the physical model. Population CMD/HB agreement, mass-function assumptions,
stellar evolution/remnants/binaries, catalogue-parameter uncertainties,
differential extinction, background estimation and Galactic potential remain
separate checks before quoting a total physical cluster mass-loss rate.

An additional physical-consistency diagnostic is
`escaped_present_day_living_mass_equivalent / adopted_current_cluster_mass =
rate * simulation_duration / adopted_current_cluster_mass` for the constant-rate
case. If this is substantial, rescaling an unchanged stream template is not a
self-consistent reconstruction: the historical progenitor mass, tidal radius
and escape kinematics may differ. This ratio mixes the declared living-star
mass equivalent with the adopted dynamical cluster mass and is a diagnostic,
not a measured historical mass fraction. The current suite does not evolve
the progenitor potential with the inferred loss rate.
