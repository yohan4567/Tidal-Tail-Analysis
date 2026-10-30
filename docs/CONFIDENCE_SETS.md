# Monte Carlo confidence sets and grid visualization

`science_validation.py` calibrates 95% **conditional finite-grid confidence
sets** using the exact template inputs and pixel grouping recorded in
`joint_grid.json`. It refuses a different template set. Observed counts,
background and signal states all use the same disjoint pixel sums.

```powershell
python science_validation.py --templates-manifest outputs/Pal5_followup/templates.json --scores outputs/Pal5_followup/scores/joint_grid.json --out outputs/Pal5_followup/validation95 --calibration 200 --validation 300 --background-scales 0.9 1.0 1.1 --seed 161803
```

The background-scale list is an explicitly declared sensitivity grid, not a
confidence interval for the background. Choose it using control-region evidence;
values between or outside those points are not covered by the stated calibration.
The command is resumable per truth/background case; changing settings, code or
inputs requires a fresh directory. Default trial counts are a starting point,
not an automatic scientific acceptance criterion.

For every candidate truth (T, rate, population scenario), simulate Poisson
counts from a randomly drawn frozen calibration state and the declared background
scale. Refit the full grid, profiling the background scale. Record separately:

- q_rate: minimum NLL at the true rate, profiling T and population;
- q_time: minimum NLL at the true T, profiling rate and population;
- q_joint: minimum NLL at true (T, rate), profiling population.

Each statistic is twice the difference from the overall minimum NLL. For m
calibration draws, use rank ceil[0.95(m+1)], or an infinite threshold if that rank
exceeds m. For each parameter coordinate, take the largest cutoff over its
declared nuisance truths. Invert these tests against the observed q values.
This avoids applying the same arbitrarily chosen asymptotic chi-square threshold
everywhere. Rate/time sets are **separate marginal sets**; their Cartesian
product is not the calibrated joint set.

The remaining held-out draws estimate coverage for every injected case; the
report includes Wilson 95% intervals for those coverage fractions. Finite-rank
coverage is marginal over repeated calibration banks and new data. It does not
guarantee 95% conditional coverage for one fixed finite bank. More draws and
independent seed repeats are needed to assess bank stability.

Outputs:

- `confidence95.json`: accepted rate and time grid values, their envelopes,
  grid-edge flags, held-out coverage and explicit limitations;
- `confidence_grid.csv`: q, calibrated threshold and acceptance for each joint node;
- `confidence_heatmap.png`: T/rate heatmap, best node and accepted joint nodes;
- `coverage.csv`: held-out coverage by injection truth and background scale.

Envelopes can contain gaps or unsampled values. They are not interpolated
endpoints. If an accepted set touches the grid edge, extend the grid before
interpreting it as a finite bound. Time plateaus need not provide a lower limit.
An empty set is reported as empty; it is not replaced with an artificial interval.

This validates frozen-template statistical inference. It does **not** validate
particle convergence, off-grid truths, a particular fixed calibration state,
misspecified population/potential, background spatial errors or absolute mass
conversion. The report deliberately leaves physical-rate validation false.
The 20+20-trial mode is only an integration/pilot test, with very uncertain coverage.

Statistical reference: [PDG Statistics review, Neyman construction and nuisance
parameters](https://pdg.lbl.gov/2023/reviews/rpp2022-rev-statistics.pdf).

The basic grid runner also saves `heatmap.png`. Cells visualize sampled nodes;
missing nodes remain blank. Irregularly spaced rates use midpoint cell edges,
not interpolation of the likelihood. White circles in the calibrated plot mark
accepted nodes rather than a visually interpolated confidence contour.

For a newly executed production grid, append `--confidence95` to
`production_study.py --execute ...` to run the matching calibration afterward.
`--calibration`, `--validation` and `--background-scales` configure that step.
It does not reuse historical fine-pixel thresholds or silently mark the physical
rate as validated.

**Compute budget:** exhaustive calibration scores the full grid for every
synthetic dataset at every truth/background point. With K candidates, B declared
background scales and m+v draws, this requires K*B*(m+v) simulated datasets and
approximately K²*B*(m+v) candidate likelihood fits. The production planner prints
this count when `--confidence95` is enabled. Start with a labeled small pilot to
check the implementation and budget. Calibrating only selected truths does not
certify coverage over an entire untested grid.
