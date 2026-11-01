# Tidal Tail Analysis

Particle-spray simulations of globular-cluster tidal tails, Gaia observation modeling, and conditional mass-loss-rate inference.

**Start with [Installation and data preparation](docs/SETUP_AND_DATA.md).** This repository contains code, tests, configuration, methods, and selected validation summaries. Gaia catalogues, the external selection project, downloaded calibration data, and full simulation outputs must be prepared separately. Paths under `outputs/` in historical experiment commands refer to locally generated products; they are not included in a fresh clone.

For fresh runs on another PC, see [NGC 6093 (M80), NGC 1851, and NGC 6341 (M92)](docs/FRESH_GC_RUN.md), including input preparation, a separate output folder, queue commands, and the NGC 1851 inference hold.

The software has been exercised on real Pal5 inputs. Its estimates remain conditional on the adopted stellar population, selection function, background, Galactic potential, and cluster evolution. A successful run or passing software test does not establish a precise physical mass-loss measurement.

Some detailed research documents are currently written in Traditional Chinese.

## Installation

Python 3.11 has been tested.

```powershell
git clone https://github.com/yohan4567/Tidal-Tail-Analysis.git
cd Tidal-Tail-Analysis
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

On Linux or macOS, use `.venv/bin/python` instead of `.venv/Scripts/python.exe`.

If `ensurepip` fails, create the environment with `python -m venv --without-pip .venv`, then use a Python installation with pip:

```powershell
python -m pip --python .venv/Scripts/python.exe install -r requirements.txt
```

The native galpy C extension could not be loaded on the development Windows machine. The tested fallback uses galpy release sampling and a Cartesian SciPy DOP853 integrator with a cached progenitor trajectory. It runs on the CPU. The native path remains available but was not validated on that machine. A C-extension warning alone does not indicate a failed simulation.

Run commands from this directory so `.galpyrc` can disable unnecessary online version checks.

## Tests without external catalogues

```powershell
New-Item -ItemType Directory -Force outputs | Out-Null
.\.venv\Scripts\python.exe -m unittest test_model test_extent test_dynamics test_mass_evolution test_joint_grid test_profile_recovery test_dense_convergence test_validation
```

On Linux or macOS, create the directory with `mkdir -p outputs`.

The public core suite contains 38 tests and also runs in GitHub Actions. The full local suite contains 60 tests, including integration checks requiring the external selection project. See [data preparation](docs/SETUP_AND_DATA.md) before attempting the complete workflow.

## External selection project

The observation pipeline reuses functions from the original `tidal_study_v2` project. Its source and catalogues are not bundled here. Point to an authorized local copy:

```powershell
$env:GC_TIDAL_STUDY = 'D:\research\tidal_study_v2'
```

```bash
export GC_TIDAL_STUDY=/path/to/tidal_study_v2
```

The legacy local layout remains a fallback:

```text
../TidaleffectsOnGCs-main/TidaleffectsOnGCs-main/tidal_study_v2/
```

Individual commands exposing `--study` also accept an explicit path. The adapter calls `step3_tail_search.py` rather than independently rewriting RUWE, flux-excess, parallax, proper-motion, CMD, neighboring-object, or missing-tile cuts. The CMD ridge and whether orbital corrections are enabled are determined from the observations, then frozen for mock selection. Outputs record source identities, checksums, and selection counts.

`--legacy` explicitly selects the original field rather than the v3 field plan. The local Pal5 example uses the available original 5-degree field; the planned 10-degree data were unavailable. There is no silent fallback. That catalogue contains 39,677 sources and yields 2,117 selected candidates, with row-by-row agreement checked against the original selector. Its narrow proper-motion download window leads to fixed proper-motion cuts; the reason is retained in the report.

## Quick start with prepared observations

After installing the dependencies and preparing the external selector and input data:

```powershell
.\.venv\Scripts\python.exe run_pipeline.py --cluster Pal5 --legacy --out outputs/Pal5_analysis
```

This default mode does not require a stream simulation or population calibration. It:

1. Reuses the original Gaia selection and masks.
2. Writes raw counts, a background map, and a predeclared two-direction orbit track.
3. Measures conditional tail-candidate segments, traceable extent, and transverse width.
4. Produces `REPORT.md`, `manifest.json`, figures, and CSV summaries.

Use a fresh output directory for each experiment. Existing results are not overwritten; failed stages retain their manifest and logs. The default mode does not report a measured physical mass-loss rate.

## Stream simulations and model scales

A small pipeline smoke test:

```powershell
.\.venv\Scripts\python.exe run_pipeline.py --cluster Pal5 --legacy --simulate --duration-myr 100 --particles 20 --out outputs/Pal5_smoke
```

Twenty particles test the interfaces only. Scientific use requires particle-density, random-seed, and integration-duration checks.

The dynamics and model-scale tools can also run separately:

```powershell
.\.venv\Scripts\python.exe stream_model.py --config configs/Pal5.json --duration-myr 1000 --particles 2000 --out outputs/Pal5_1gyr
.\.venv\Scripts\python.exe model_extent.py --run outputs/Pal5_1gyr --out outputs/Pal5_model_extent
```

Shared assumptions include:

- galpy 1.12 `chen24spraydf`, with equal numbers of leading- and trailing-arm tracers.
- `MWPotential2014`, normalized to `ro = 8 kpc` and `vo = 220 km/s`.
- A moving Plummer cluster potential with `a = r_hm * sqrt(2^(2/3) - 1)`.
- galpy sky transformations using `z_sun = 20.8 pc` and solar peculiar velocity `[-11.1, 12.24, 7.25] km/s` in its convention.
- Tracer mass weights representing both arms together, without an extra factor of two.
- No tracer self-gravity or explicit internal relaxation. The external selector retains its own coordinate and orbit assumptions, which require sensitivity checks.

### Fixed-dynamics baseline

The supplied `configs/Pal5.json` uses a fixed cluster mass and uniform release-time sampling. Each tracer represents `reference_mdot * duration / n_particles`. Changing the rate amplitude changes weights rather than the dynamics. This approximation needs reassessment when the represented escaped mass is substantial.

### Prescribed evolving-mass models

An explicit `mass_evolution` configuration enables a prescribed mass history used consistently in release kinematics and post-release cluster forces. The joint-study driver uses stratified release times and can reuse a maximum-duration simulation for shorter windows under the documented conditions. Rates that change the potential require separate dynamical models; a free amplitude cannot substitute for those models.

This is a prescribed history, not a self-consistent prediction of relaxation, stellar evolution, radius evolution, or the escape rate. See [Joint time and mass-loss experiments](docs/JOINT_GRID.md) and the [experimental method](docs/EXPERIMENTAL_METHOD_ZH.md).

Simulation products include `particles.csv`, `orbit.csv`, `overview.png`, and `run.json`; applicable modes also save mass histories, released states, and integration checkpoints. Particle columns describe present-day positions, distances, proper motions, radial velocities, release lookback times, mass weights, and model arm labels. Negative orbit times indicate the past; Cartesian axes follow galpy conventions.

`model_extent.json` reports per-arm 95%-mass angular radial scales and model 3D distances from the cluster. These are quantiles, not along-tail lengths, observed detection limits, or Jacobi radii. Fewer than 100 particles per arm triggers a small-sample flag; exceeding that count does not establish convergence.

## Observed tail extent and width

```powershell
.\.venv\Scripts\python.exe tail_extent.py --selection outputs/Pal5_analysis/selection --out outputs/Pal5_profile --bin-deg 0.4 --half-width-deg 0.3
```

The method:

- Defines an orbit from catalogue phase-space coordinates, rather than choosing the brightest direction in the same density map.
- Counts unsmoothed observations and the original background in disjoint bins along the projected track.
- Excludes the central `1.5*r_t` by default and searches to `0.9*min(field radius, 5 deg)`.
- Uses a default bin width of `max(0.25 deg, 3 pixels)` and transverse half-width of `max(0.25 deg, 2.1 pixels)`, both configurable.
- Applies fixed-background Poisson upper-tail tests with Bonferroni correction over valid bins in both directions.
- Requires at least 70% usable strip area and two consecutive candidate bins by default. Masks and track or field boundaries interrupt the connected extent.
- Reports connected extent and separate candidate segments. No signal yields a null extent rather than a zero-length measurement.

The extent bracket expresses bin resolution, not a confidence interval. Width is a background-subtracted transverse second-moment RMS within candidate bins, affected by the strip boundary and pixels; it is not a deconvolved Gaussian dispersion.

Forward/backward labels follow the projected orbit and do not automatically identify physical leading/trailing arms. Angular distances are projected arc lengths. Conversion to parsecs uses the cluster distance and is not a 3D deprojection. A missing connected segment does not establish the absence of a tail; this is not a general halo, shell, or off-orbit structure detector.

The background is estimated from the same field and may absorb signal. Background-shape uncertainty is not fully propagated. Candidate p-values are conditional on the fixed background and orbit and do not cover searches over multiple directions, widths, settings, or clusters. Exploratory `sqrt(N+B)` error bars are not a complete background uncertainty model.

## Injection recovery and false signals

```powershell
.\.venv\Scripts\python.exe validate_tail.py --selection outputs/Pal5_analysis/selection --true-extent-deg 2 --selected-stars 300 --sigma-deg 0.1 --trials 100 --out outputs/Pal5_recovery
```

Here, 300 is the total expected selected count in both arms within usable sky, not a mass or rate. This experiment draws Poisson background-only and injected catalogues against a fixed background map and reports detection fractions, binomial intervals, recovered endpoints, and bias. It does not rebuild the CMD ridge, refit the background, or validate survey completeness.

A historical fixed-seed local test found no candidates in 100 background-only trials, with an approximately 3.6% upper 95% binomial bound. All 100 injections with 300 selected stars and a 2-degree extent produced connected candidates, with endpoint bias roughly within one bin. These are conditional software experiments, not guarantees of the real-sky false-alarm rate or completeness. Test a range of strengths, widths, lengths, and settings.

## Population and Gaia observation modeling

The recommended observation path is **Hybrid v2**: official TCG completeness maps, population quadrature, availability and quality probabilities, correlated astrometric errors, weighted templates, and shared calibration-state sampling.

- [Observation workflow](docs/OBSERVATION_V2.md)
- [Population and empirical Gaia calibration](docs/SCIENCE_CALIBRATION.md)
- [Validation workflow](docs/VALIDATION.md)

`examples/create_v2_demo.py` and the small `*_SYNTHETIC.csv` files are software fixtures. They must not be used as empirical calibration or to claim an observed physical mass-loss rate.

### Legacy calibrated mode

The older `--mode calibrated` path remains as a donor-based comparison. It requires:

| Input | Meaning |
|---|---|
| `population.csv` | Absolute G, intrinsic BP-RP, and stars per unit represented escaped mass |
| `completeness.csv` | Apparent G and detection/catalogue-entry probability before scientific quality cuts |

Define the population mass basis, including the treatment of unseen stars and remnants. Completeness must cover the simulated magnitude range and must not double-count later quality selection. Supply extinction, reddening, and their sources explicitly.

```powershell
.\.venv\Scripts\python.exe run_pipeline.py --cluster Pal5 --legacy --mode calibrated --reuse-run outputs/Pal5_smoke/simulation --population population.csv --completeness completeness.csv --a-g 0.1 --reddening 0.05 --calibration-description "Describe population and completeness sources" --null-trials 200 --out outputs/Pal5_calibrated
```

The extinction values illustrate syntax only. Replace them with justified inputs and replace the smoke simulation with an adequately sampled model. Synthetic fixtures require `--synthetic-calibration` and remain explicitly labeled.

This legacy approximation samples luminosity bins, applies completeness weights, draws error and quality donors from the original field, perturbs astrometry, and reuses the original cuts. It does not provide the full covariance, photometric-scatter, or spatial-systematics treatment of Hybrid v2. More tracer rows do not automatically mean more physical stars.

Its model is `N_j ~ Poisson(A*S_j + beta*B_j)`, with nonnegative signal and background amplitudes. The reported rate is `A*reference_rate`. Asymptotic profile intervals and upper limits are approximate near boundaries or at low counts and do not include all population, completeness, or template errors. A fixed-template null bootstrap is not a blind-search significance calculation. Unidentifiable templates or missing signal are reported as failures rather than manufactured rates.

A standalone synthetic count example:

```powershell
.\.venv\Scripts\python.exe fit_rate.py --demo --null-trials 100 --out outputs/synthetic_rate_demo.json
```

## Joint time and mass-loss inference

The coupled constant-rate workflow supports stratified sampling, maximum-time window reuse, background profiling, discrete population scenarios, independent-template diagnostics, and frozen-template recovery. It does not automatically produce a calibrated continuous confidence interval or a validated total physical rate.

| Document | Scope |
|---|---|
| [Joint-grid guide](docs/JOINT_GRID.md) | Coupled runs, window reuse, scoring, and recovery |
| [Experimental method](docs/EXPERIMENTAL_METHOD_ZH.md) | Definitions, assumptions, and reporting levels |
| [Initial joint validation](docs/JOINT_VALIDATION_RESULTS.md) | Passed and failed numerical checks |
| [Improvement research](docs/IMPROVEMENT_RESEARCH_ZH.md) | Literature review and resolution diagnostics |
| [Dense sampling and profile tests](docs/DENSE_TEST_RESULTS.md) | Three-seed 8,000-particle runs, a sparse 4-Gyr preflight, and enlarged independent calibration |
| [Selected aggregate results](docs/results/README.md) | Small public tables and figures |

The dense tests improved numerical stability, but fine spatial templates did not pass every target. Profile-statistic calibration with exact frozen templates does not validate population assumptions, Galactic dynamics, background construction, or the complete selection pipeline.

## Main tools

| Tool | Purpose |
|---|---|
| `run_pipeline.py` | Pipeline stages, logs, manifests, and reports |
| `stream_model.py` | Chen particle release and orbital integration |
| `make_config.py` | Import a cluster from an external catalogue |
| `selection_bridge.py` | Reuse observation and mock selection |
| `response_model.py` | Hybrid population and Gaia forward model |
| `run_joint_study.py` | Coupled rate grid with reusable time windows |
| `joint_grid.py` / `joint_fit.py` | Conditional scores and globally bounded background profiling |
| `profile_recovery.py` | Separate rate/time profile calibration with explicit seeds |
| `dense_convergence.py` | Independent-seed and nested-density diagnostics |
| `tail_extent.py` / `model_extent.py` | Observed candidate extent and model quantile scales |
| `validate_tail.py` | Selected-space injection and null tests |
| `fit_rate.py` | Fixed-dynamics amplitude inference and conditional intervals |

Full integration tests with the external selector available:

```powershell
.\.venv\Scripts\python.exe -m unittest test_profile_recovery test_dense_convergence test_joint_grid test_mass_evolution test_validation test_science_calibration test_observation test_model test_selection test_extent test_dynamics
```

They cover normalization, weights, selection freezing, masks, projection, extent recovery, background-only cases, profile statistics, finite calibration ranks, checkpointing, and dynamical consistency. Observation selection is additionally checked against the original selector during each selection run.

Scientific use still requires justified population and completeness inputs, background and spatial-error checks, parameter and seed sensitivity, end-to-end injection recovery, and Galactic-potential and cluster-evolution comparisons. See the [implementation plan](docs/IMPLEMENTATION_PLAN.md) for broader research scope; current behavior is documented in the workflow guides.

## References

- Chen et al., [Improved particle spray algorithm for modeling globular cluster streams](https://arxiv.org/abs/2408.01496).
- [The authors' particle-spray examples](https://github.com/ybillchen/particle_spray).
- [galpy documentation](https://docs.galpy.org/).

The implementation calls the relevant libraries rather than copying the authors' particle-spray code. Further methodological references are provided in the research documents.

## License

This project's code is available under the [MIT License](LICENSE).
External dependencies, separately obtained selection code, and datasets retain their respective licenses and terms.

## Higher-resolution follow-up grid

`production_study.py` supplies a reproducible follow-up preset: 16 rates
(1 to 31 solar masses/Myr, step 2), 20 lookback windows
from 100 to 2000 Myr (step 100), three independent dynamics seeds, and a 0.25 Myr release
interval. Each arm receives one tracer per interval. Each Tmax simulation has
16,000 tracers; the full preset integrates 768,000 trajectories across 320 T/rate nodes. This is a
substantial computation; no runtime guarantee is implied.

```powershell
python production_study.py --config PATH_TO_COUPLED_DYNAMICS_JSON --selection PATH_TO_SELECTION --response-configs POP_MINUS1_JSON POP_ZERO_JSON POP_PLUS1_JSON --cluster-age-myr 11000 --out outputs/Pal5_followup
```

The command prints a plan. Append `--execute` to execute the resumable grid.
Supply the actual externally justified age bound and calibrated population
configurations; the example age is not an age measurement.

All candidates use the same particle density and disjoint 4-by-4 pixel sums.
Observed integer counts, background expectations, and every signal calibration
state are summed with the same map. Masked pixels remain excluded. This is
coarser spatial binning, not smoothing observed counts. The factor was motivated
by the separate Pal 5 density experiment; other clusters need their own check.
`run_joint_study.py --pixel-factor` also exposes this setting explicitly.

The assessment reports boundary optima and outstanding validation requirements.
It never certifies an absolute physical rate merely because this preset ran.
The grid remains conditional on population, mass conversion, potential and
selection assumptions. Extend the rate grid if its optimum touches an edge;
interpret a flat time profile as loss of identifiability, not an automatic lower
limit. Independent refinement and recovery tests, including weak signals, grid
edges and short/long windows, remain necessary.

**Calibration compatibility:** existing fine-pixel recovery thresholds must not
be applied to the aggregated likelihood. The current `profile_recovery.py`
consumes original templates; it does not automatically reproduce the production
pixel aggregation. Its historical thresholds therefore do not certify this
preset. No confidence interval is emitted by the production assessment until a
matching recovery analysis is available. Earlier pilot outputs remain pilot
outputs and are not overwritten.

## Monte Carlo 95% confidence sets

The matching-bin validation runner and heatmap output are documented in
[Confidence sets and visualization](docs/CONFIDENCE_SETS.md). It calibrates rate,
time and joint T/rate sets separately, saves held-out coverage diagnostics and
marks grid-edge acceptance. These remain conditional finite-grid statistical
results; they do not certify the physical population model or absolute rate.
This supersedes the earlier limitation that production aggregation had no
matching recovery implementation. Historical fine-pixel calibration thresholds
still must not be reused.


The grid axes are configurable independently, for example:

```powershell
python production_study.py --config PATH_TO_COUPLED_DYNAMICS_JSON --selection PATH_TO_SELECTION --response-configs POP_ZERO_JSON --cluster-age-myr 11000 --out outputs/Pal5_dense_grid --rate-min 1 --rate-max 31 --rate-step 2 --time-min-myr 100 --time-max-myr 2000 --time-step-myr 100
```

This prints the plan; append `--execute` to run it. Smaller time steps reuse
Tmax trajectories; additional rate nodes require independent dynamics. Dense
sampling does not repair missing support beyond an accepted grid edge. Expand
bounds separately if needed. No interpolation of the old 3-by-4 pilot is used.
Exhaustive per-node Monte Carlo calibration grows approximately quadratically
with the number of candidates; inspect the printed calibration cost first.

## Accelerated dynamics and restart

The batched backend, numerical validation, 10-worker Pal 5 command and
checkpoint/restart behavior are described in [Performance](docs/PERFORMANCE.md).
The accelerated run uses stricter integration tolerances and keeps the old
reference outputs separate. Scientific validity is still conditional on the
population, selection and dynamical assumptions.


When changing a grid, `--reuse-grid PATH_TO_INACTIVE_GRID` can copy matching
completed simulations and integration checkpoints. Reuse requires identical
physical/numerical configurations, integration source hashes, package versions
and file checksums. A checkpoint's complete input digest is checked again by the
integrator on resume. Keep the source grid inactive and retain its original
files. Incompatible simulations are refused rather than silently reused.
The Pal 5 launcher now selects 1, 3, ..., 31 solar masses/Myr and uses
`outputs/cluster_survey_29_v1/Pal5/production_grid_shared_v1`; the inactive odd-rate
grid is its verified dynamics reuse source. Compatible population and quality
responses now share observation draws by default while retaining separate outputs
and all calibration samples. Backend provenance and cache keys distinguish the
new implementation. See [Performance](docs/PERFORMANCE.md) for equality checks and
limits on interpreting the measured speedup.


## Continue the cluster sample

The resumable remaining-cluster queue, cluster-specific population assumptions,
input checks and progress files are documented in [Survey execution](docs/SURVEY_EXECUTION.md).
The completed Pal 5 grid is conditional; see the [literature cross-check](docs/results/Pal5_literature_crosscheck.md)
before interpreting its boundary optimum as a physical mass-loss measurement.


## Population-dependent escape rates

Production reports now separate the conditional rate profiles for each population
scenario and leave the primary absolute rate unset while the escaping population
normalization lacks external calibration. [Population diagnostics](docs/POPULATION_DIAGNOSTICS.md)
adds visibility-conversion curves and spatial residual maps using existing scored
templates; likelihood weights and previously computed scores are unchanged.


## Marginalized posterior and measured-MF populations

Production post-processing now integrates the raw Poisson-mixture likelihood over
background and time, reporting conditional posterior densities and credible
intervals with explicit prior sensitivity. Legacy population scenarios remain
separate unless a justified population prior is supplied. The empirical pipeline
uses measured global mass-function slopes and separate escaping-MF transfer
assumptions. See [Bayesian population inference](docs/BAYESIAN_POPULATION_INFERENCE.md)
for priors, phase weights, numerical tests, running jobs and scientific limits.

## Directly measured tail populations

`tail_population.py` can construct a population quadrature from a published tail
mass-function measurement instead of assuming an escape offset from the global
cluster MF. `configs/Pal5_tail_mf_Koch2004.json` records the measured north/south
segments separately, including the slope convention and fitted mass interval.
The builder refuses extrapolation beyond that interval by default. An explicit
`--allow-extrapolation` creates a **conditional** model, never a validated total
absolute escape rate. It does not overwrite existing grids or start simulations.
See [absolute-rate requirements](docs/ABSOLUTE_RATE_REQUIREMENTS_ZH.md) for the
remaining mass-budget, selection and dynamical constraints. More precise
conditional credible intervals cannot replace missing physical information.

The direct empirical runner omits the redundant legacy tilt-response grid while
keeping the existing rates, time windows, seeds and tracer counts. A serial
cluster queue uses up to twelve workers per cluster and records input blockers;
see [performance and queue details](docs/PERFORMANCE.md#direct-empirical-queue-2026-10-31).

The [workflow audit](docs/EXPERIMENT_AUDIT_ZH.md) distinguishes conditional rate
inference from total mass-budget closure. New posteriors include a background-only
model comparison with properly normalized finite axis priors; this is not a
calibrated detection significance. `run_tail_study.py --workers 8` reruns Pal5
using separate north/south literature MF transfer hypotheses. Other clusters'
population constraints and remaining data needs are described in the
[escaping-population strategy](docs/ESCAPING_POPULATION_STRATEGY_ZH.md).

All 29 targets now have a sourced [population registry](docs/POPULATION_REGISTRY_ZH.md).
`configs/population_supplements_v1.json` supplies additional literature ages while
preserving frozen inputs and multi-population holds. The queue loads these
supplements. `python prepare_population_inputs.py` downloads only missing official
isochrones and records failures without overwriting existing files; `python
population_catalog.py` refreshes the registry and its model-validity checks.
An indirect global MF estimate (including ESO452-SC11's catalogue value) is not
accepted as a directly measured escaping or global MF constraint.

## Continuous survey and execution logs

On Windows, run `powershell -File scripts/run_continuous_survey.ps1` after preparing
the external data. The default is one cluster at a time with 12 numerical/response
workers, single-threaded BLAS/OpenMP per worker, and a 24 GiB soft RAM budget.
Successful, failed, or blocked clusters advance to the next target automatically.
Restarting prioritizes interrupted work and resumes compatible checkpoints;
completed results are skipped. Pal 5's cancelled rerun remains excluded.
Known failures are not repeatedly retried; use `-RetryFailed` after correcting them.

The RAM setting checks system headroom before starting a stage; it is not an OS
allocation limit and does not allocate unused memory. On other machines, adjust
`-MemoryBudgetGiB` and `-Workers`. Scientific grids, seeds, integration batches,
and frozen provenance are unchanged. Background quadrature is still serial.

Read `outputs/cluster_survey_29_v1/EXECUTION_LOG.md` for the current pass,
`execution_events.jsonl` for append-only commands, timestamps, exit codes and
30-second system-RAM samples, and `PROGRESS.md` for numerical progress. Per-cluster
logs retain full errors. All targets being attempted does not mean they succeeded
or passed scientific validation. Do not run two writers against the same root.
