# Assumed escape rate versus relative Poisson fit

The calibrated pipeline automatically writes `rate_profile/profile.png`, an SVG,
the scanned values in CSV, and metadata in JSON. To plot an existing response:

```powershell
.\.venv\Scripts\python.exe plot_rate_profile.py --input outputs/Pal5_hybrid_v2/response/calibrated_counts.npz --out outputs/Pal5_hybrid_v2/rate_profile_new --label Pal5
```

Use a new output directory. Optional `--points`, `--min-rate`, and `--max-rate`
control the scan; the range must contain the fitted optimum.

The horizontal axis is the assumed escaped stellar mass per unit time. The
vertical axis is q = -2 [ln L(rate) - ln L(best)], so lower is better. Each rate
refits the background normalization with its shape fixed. Shared response
calibration draws are averaged at the whole-data likelihood level. For fixed
dynamics and escape history shape, rates rescale the existing signal template;
they do not require rerunning the particle trajectories.

The shaded region uses the local asymptotic q = 1 interval. The q = 2.71 line
marks the approximate one-sided 95% upper-limit threshold on the upper branch.
These approximations need calibration near boundaries or for multimodal curves.
This curve compares assumed rates within the chosen model; it is neither an
absolute goodness-of-fit test nor a p-value or posterior distribution.

The supplied Pal5 hybrid demonstration has synthetic population and response
calibration. Its numerical optimum is a software demonstration, not a measured
physical Pal5 mass-loss rate. Unknown selection, population, background-shape,
and dynamical-model systematics are not covered by the plotted intervals.
