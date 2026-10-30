# Population-conditioned reporting and spatial checks

The grid's MLR normalization is conditional on its escaping-population model.
A global optimum over uncalibrated population variants is retained as a numerical
candidate, but is not reported as the primary absolute rate. Production assessments
now explicitly set `primary_rate_msun_myr: null` and include conditional optima,
boundary flags and their sensitivity span. This span is not a confidence interval
or a calibrated systematic-error bound. It may underestimate sensitivity when an
optimum reaches the grid edge.

`population_diagnostics.py` is post-processing. It leaves the original grid scores,
likelihood weights, dynamics and response templates unchanged. It reconstructs
selected predictions from the exact scored input hashes, uses the same disjoint
pixel aggregation, and averages independent dynamics realizations with aligned
calibration states. Native arrays and a frozen selection-map hash are checked.

Artifacts in `<grid>/population_diagnostics_v1`:

- `scenario_best.csv`: fixed-scenario optima and global relative q.
- `conditional_rate_profiles.csv` / `population_profiles.png`: profile over time,
  separately re-zeroed within each full population/calibration scenario. Curves
  cannot be compared vertically as global likelihood differences; use the table.
- `population_visibility_grid.csv` / `population_visibility.png`: expected selected
  counts and selected stars per represented escaped mass at all grid nodes. The
  plot uses a common window closest to half the maximum simulated time, chosen
  without reference to fit quality.
- `candidate_*.npz` / residual PNGs: observed counts, fitted background, signal,
  total prediction and descriptive Pearson residuals for each scenario optimum
  and the best available candidates at rates 5 and 31. Those comparison rates are
  diagnostic examples, not priors or special constraints in the fitting.
- `population_summary.json` / `REPORT.md`: provenance, definitions and limits.

Visibility is N_selected / (rate * T) for the modeled window. It includes the
stellar population, dynamics, field, lookback interval and selection. It is not a
universal stars-per-mass correction and cannot be held fixed while changing the
rate/time without checking the model. The expected escaped mass is the represented
living-star mass, not a complete accounting of dark remnants or stellar winds.

Spatial residuals use the calibration-mean prediction, not a posterior-weighted
calibration prediction. Their mean-prediction Poisson deviance/Pearson sums are
named explicitly and are not the mixture-likelihood q, calibrated significance,
or an automatic goodness-of-fit test. No p-value or confidence boundary is drawn.

Run on completed data:

```powershell
python population_diagnostics.py --grid PATH_TO_GRID --selection PATH_TO_SELECTION
```

Future production runs invoke this after scoring. A separate completion watcher
can process already running batches without changing their source-frozen scientific
code: `python population_diagnostics.py --survey-root outputs/cluster_survey_29_v1`.
It waits for the production assessment and checks every 30 seconds by default.
The output directory is reserved atomically to avoid two diagnostic writers; a
summary completion marker is written last. Incomplete directories/errors require
inspection. Existing assessments are preserved inside the new diagnostic folder
before their reporting fields are updated; original best-candidate values remain.

## What this change does not establish

No external luminosity/mass-function constraint has yet been applied. It does not
measure a new physical MLR, force Pal 5 toward a literature value, change q weights,
or break a physical normalization degeneracy through plotting alone. Constraining
the escaping luminosity function requires a background-aware stellar-population
likelihood and a validated observational selection model. Reusing the same selected
stars in an extra likelihood term must not count their evidence twice. A core mass
function needs a justified dynamical transfer model before being treated as the
escaping mass function; unseen low-mass stars remain dependent on extrapolation.

For Pal 5 the three conditional grid optima are 31, 11 and 5 solar masses/Myr for
tilts -1, 0 and +1. The global and rate-5 candidates predict about 85 and 100 selected
stars in their own best windows, despite very different represented escaped masses.
These differences motivate external normalization and background checks, not an
ad hoc preference for either answer. See the separate literature cross-check.
