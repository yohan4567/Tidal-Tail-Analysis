# Real population and empirical Gaia observation calibration

This workflow replaces the integration-test fixtures. It does not label a model
as scientifically validated merely because its inputs are no longer synthetic.
The mass-loss inference remains conditional on population, extinction, Galactic
potential, background and selection assumptions.

## Sources and mass definition

* Dynamics: the existing project's `catalogs/master_clusters.csv` snapshot.
  `make_config.py` records its hash and imported columns. Positions/velocities
  originate from Baumgardt's `orbits_table.txt`; mass and three-dimensional
  half-mass radius originate from `combined_table.txt`.
* Pal5 population baseline: age 11 Gyr and scaled-solar [M/H] = -1.4. This is a
  conditional choice motivated by the 11 Gyr, -1.4 population used by Balbinot &
  Gieles (2018), not a new determination of Pal5's age or chemical composition.
  [M/H] must not silently be equated to spectroscopic [Fe/H] for alpha-enhanced
  populations. The age, abundance and photometric CMD agreement need sensitivity
  checks before a scientific result is quoted.
* Official PARSEC v1.2S/COLIBRI, CMD 3.8, Gaia EDR3 passbands and its Kroupa
  cumulative IMF. Requests, returned file hashes and model assumptions are saved.
  Downloaded stellar-evolution predictions are not arbitrary synthetic fixtures.
* The baseline is renormalized to the **living stellar mass represented by the
  supplied isochrone**. Artificial dark termination/remnant points (PARSEC label
  9 or logL <= -9) are excluded. This is not automatically the total cluster
  dynamical mass, total initial mass, or a complete remnant mass census.
* Escaping stars need not follow the initial mass function. Three explicit
  sensitivity populations weight the baseline by `(current_mass / 0.5)**tilt`,
  tilt = -1, 0, +1, then renormalize the living stellar mass to one solar mass.
  They preserve the isochrone evolutionary-phase weights. These are hypotheses,
  not a fitted escape mass function, pure power-law slopes, or confidence bounds.
  A mass-function change can be nearly degenerate with the rate amplitude;
  similar likelihood values do not establish which escape mass function is true.

## Commands

Run from `gc_stream_model` in PowerShell; use fresh output paths. Outputs and
downloads are under the ignored `outputs/` directory.

```powershell
.\.venv\Scripts\python.exe fetch_isochrone.py --age-years 11000000000 --mh -1.4 --source "Pal5 conditional baseline; Balbinot and Gieles 2018 MNRAS 474 2479; scaled-solar abundance assumption" --out outputs/Pal5_science_calibration/isochrone.dat

.\.venv\Scripts\python.exe selection_model.py parent-query --ra 229.02 --dec -0.121 --radius 5.05 --g-min 12 --g-max 21 --sampling-modulus 10 --download --max-records 150000 --out outputs/Pal5_science_calibration/parent_random10.csv.gz

.\.venv\Scripts\python.exe prepare_science_calibration.py --isochrone outputs/Pal5_science_calibration/isochrone.dat --parent outputs/Pal5_science_calibration/parent_random10.csv.gz --selection outputs/Pal5_hybrid_v2/selection --out outputs/Pal5_science_calibration/bundle

.\.venv\Scripts\python.exe run_pipeline.py --cluster Pal5 --legacy --mode weighted --reuse-run outputs/Pal5_hybrid_v2/simulation --response-config outputs/Pal5_science_calibration/bundle/tilt_plus0_zp_minus0.017.json --out outputs/Pal5_real_inputs_baseline
```

The downloader supports `--endpoint aip` for the AIP Gaia DR3 mirror. It rejects
TAP ERROR/OVERFLOW and incomplete files, and never treats the first N sources of
a truncated download as a representative parent. `MOD(random_index,10)=0` is an
independent random subsample spanning the whole sky footprint, not a PM/CMD or
quality cut. It estimates conditional quality probabilities and error donors;
do not use its raw size as the total field-star count. A representative subset
can be used for these conditional probabilities without multiplying the
estimated probabilities by the inverse sampling fraction.

If cone queries or transfers fail, use the resumable indexed source-ID route;
each disjoint patch is separately checked for TAP overflow and hashed:

```powershell
.\.venv\Scripts\python.exe fetch_parent_tiles.py --ra 229.02 --dec -0.121 --radius 5.05 --sampling-modulus 10 --out outputs/Pal5_science_calibration/parent_tiled.csv.gz
```

Pass this parent path to `prepare_science_calibration.py`. Successful patch
downloads are reused on the next invocation. An interrupted patch is not data.

## What is empirically calibrated

The parent includes raw G/BP/RP photometry, RUWE, flux excess, astrometric errors,
their reported correlations and flux SNRs. Null entries are retained in the
denominator. No membership or quality cuts are applied at download time.
The original PM-selected project field remains the observational target sample;
it must not be substituted for this parent when learning quality probabilities.
Preparation excludes the cluster core and masked neighbouring clusters by sky
position alone, since those crowded regions are excluded from the science fit.
Both successes and denominators get the same exclusions, recorded in the donor
parent metadata. This reduces contamination of the tail-field calibration by
crowded cluster cores; it does not select by PM/CMD/quality.

The code uses a fixed sky/G/G-RP grid, splits 80/20 before learning the map, and
saves held-out availability/quality Brier scores, predictive residuals and
coverage. It then refits the map with the full parent. Error vectors are donated
jointly by nearby real sources in magnitude, colour and sky position; the
reported astrometric covariance is preserved. It does not draw arbitrary RUWE
values or pretend all faint stars have the same error. Unsupported calibration
cells or missing local error donors raise errors. Held-out field-star tests do
not by themselves prove transferability to the escaping population.

## Remaining declared nuisance assumptions

The supplied bundle has nine scenarios: three population tilts times scalar
parallax offsets -0.030, -0.017 and 0 mas. These offsets are **sensitivity
settings**, not a full per-source Gaia zero-point correction or uncertainty
distribution. The baseline E(B-V)=0.03 comes from the project's Harris2010
snapshot. Band extinctions use CMD's G2V coefficients (G/BP/RP =
0.83627/1.08337/0.63439 times Av, with Av=3.1 E(B-V)); these are a low-extinction
constant-field approximation, not a colour-dependent or differential dust map.
Scan extinction, age, composition and mass limits as additional models when
reporting systematics. Scenario spread must not be labelled a statistical CI.

Before reporting a physical rate: inspect the population CMD against independent
cluster members; review held-out calibration and sky/colour support; establish
particle, photometric quadrature and Monte Carlo convergence; check alternate
backgrounds and dynamics; and calibrate goodness-of-fit with mocks. The old
`Pal5_hybrid_v2` outputs retain their SYNTHETIC labels and are not overwritten.

The current CMD diagnostic uses quality/astrometry-selected core stars **before**
the empirical CMD cut. The initial fixed-Reimers PARSEC baseline does not reproduce
Pal5's blue horizontal branch. It therefore needs an HB prescription/sensitivity
check or a consistently restricted science sample before a physical rate is
adopted; simply having an official model file does not resolve this mismatch.
The conversion also assumes the same present-day population at all release
times and no mass-dependent escape geometry. Its normalization is present-day
living-stellar-mass equivalent, not a full stellar-evolution mass budget at each
historical escape epoch.

```powershell
.\.venv\Scripts\python.exe check_population_cmd.py --population outputs/Pal5_science_calibration/population_0 --selection outputs/Pal5_hybrid_v2/selection --out outputs/Pal5_science_calibration/cmd_check_new
```

References:

* PARSEC CMD: https://stev.oapd.inaf.it/cgi-bin/cmd_3.8
* Population and stream evolution: https://academic.oup.com/mnras/article/474/2/2479/4562632
* Gaia subsample selection: https://doi.org/10.1051/0004-6361/202346547
* Gaia DR3 survey completeness: https://arxiv.org/abs/2208.09335
