# Tidal Tail Analysis

Particle-spray simulations of globular-cluster tidal tails, Gaia observation modeling, and conditional mass-loss-rate inference.

**Start with [Installation and data preparation](docs/SETUP_AND_DATA.md).** This repository contains code, tests, configuration, methods, and selected validation summaries. Gaia catalogues, the external selection project, downloaded calibration data, and full simulation outputs must be prepared separately. Paths under `outputs/` in historical experiment commands refer to locally generated products; they are not included in a fresh clone.

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
