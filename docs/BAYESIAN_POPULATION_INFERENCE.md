# Marginalized inference and empirical population models

The new inference reconstructs raw Poisson-mixture predictions from scored,
checksummed templates. For every rate/time/population node it integrates the
background normalization rather than substituting its maximum-likelihood value:

    p(rate | data) proportional to sum_k w_k integral L(data | rate,T,beta,k)
        p(rate) p(T) p(beta) dT d(beta)

Calibration states already represent an empirical Beta-cell calibration mixture;
independent dynamics realizations are averaged before evaluating this mixture.
The code does not exponentiate a profiled q curve and label it a marginalized
posterior. `fixed_rate_fit` supplies a numerical scaling constant only: the actual
integration still evaluates the raw Poisson likelihood over beta.

## Explicit priors and outputs

- Background beta: uniform on [0,5], chosen independently of the fitted MLR.
  `--beta-upper` permits separate prior-bound sensitivity runs.
- Rate/time: the finite support of the supplied grid. Uniform and log-uniform
  alternatives are reported in all four combinations, with nonuniform-grid
  trapezoidal integration weights. No extrapolation beyond the grid occurs.
- Without a justified population prior, legacy tilt scenarios get separate
  conditional posteriors. They are not mixed with arbitrary equal probabilities.
- With the empirical population prior, Gaussian-quadrature weights integrate
  measured global-MF uncertainty. These are quadrature weights, not hand-tuned
  q weights or equal model probabilities.

`bayesian_grid.py` saves integrated node likelihoods, rate/time joint and marginal
densities, posterior medians, central 68% and 95% credible intervals, sampled mode
locations, endpoint quadrature masses and prior-sensitivity plots. Quantiles are
computed from a piecewise-linear density, not a discrete-node CDF. Endpoint masses
are resolution/support diagnostics; they are not unobserved probability beyond
the grid. Wide/multimodal distributions and upper-edge modes must not be hidden by
reporting a median alone. These intervals are Bayesian credible intervals, not
frequentist confidence intervals or a validation of physical accuracy.

Node integrals are checkpointed and the posterior output directory freezes its
input/prior/code plan. A different plan requires a fresh directory. Inspect an
incomplete checkpoint line after a crash rather than silently accepting it.
A background-only Poisson example is tested against its analytic truncated-Gamma
integral, and grid priors/quantiles have separate numerical tests.

## Empirical mass-function construction

`empirical_population.py` uses a saved official Baumgardt catalogue HTML snapshot,
records its hash, extracts fitted mass range, slope/error and checks the green
indirect-estimate flag. Indirectly estimated MF slopes are not silently promoted
to direct observations. The convention is **dN/dm proportional to m^a** (Salpeter
has a=-2.3), rather than dN/dm proportional to m^(-a).

Verified inputs:

| Cluster | Global MF slope | Reported error | Measured mass range (Msun) |
|---|---:|---:|---|
| Pal 5 | -0.84 | 0.19 | 0.29--0.81 |
| NGC288 | -0.66 | 0.11 | 0.22--0.82 |

The measured global MF is used as a Gaussian conditional input distribution.
The baseline assumes the escaping MF has the same slope. Separate offsets of
0, -0.5 and -1 explore progressively greater relative escape of low-mass stars.
**Those offsets are sensitivity assumptions, not measured probabilities**; their
posteriors remain separate. Five Gauss-Hermite nodes represent the measured slope
uncertainty. Increase node order with a fresh `--study-label` (for example, order7_v1) to assess population-quadrature
convergence; the initial five-node run is not itself evidence of convergence.

The population uses a monotone main-sequence initial/current-mass map and
power-law number integrals in unevolved stellar mass. Evolved phases retain their
narrow progenitor-mass intervals along the isochrone; giant rows are not assigned
equal probabilities and their wind-reduced present mass does not become a new
birth mass. Every population is renormalized to one solar mass of its represented
living stars, with unit-mass tests and explicit source/phase assumptions.

The single slope is extrapolated outside the measured mass range. Photometry is
not invented below the isochrone's lower mass limit; that omitted low-mass range,
remnants, binaries and winds are not magically recovered. The result remains a
rate for the explicitly represented stellar mass. Age, metallicity/alpha-element,
extinction, orbit, potential and background-shape systematics remain conditional.
Global cluster MF data are an external input, not a measured escaping MF.

## Execution

Existing templates can be reanalysed without redoing dynamics:

```powershell
python bayesian_grid.py --grid PATH_TO_GRID --out NEW_POSTERIOR_DIRECTORY
```

New production runs now invoke this inference after their existing diagnostics.
The empirical pipeline saves the official MF HTML snapshot if it is absent, then
prepares cluster-specific constrained populations, reuses
verified 48-run dynamics, rebuilds response templates and evaluates posteriors:

```powershell
python run_empirical_study.py --clusters Pal5 NGC288 --response-workers 2
```

Outputs are separate under `<cluster>/empirical_mf_v1/offset_*/`; the live state is
`<cluster>/empirical_job.json`, with global log `empirical_study.stdout.log`.
The old arbitrary-tilt grids and outputs are retained for comparison. Completed
Pal 5 and NGC288 source grids are the first supported inputs; other clusters
need their own verified calibration, literature MF and dynamics. No partial or
failed cluster is substituted by Pal 5 data. The new posterior job must finish
before any new measured-MF MLR summary can be quoted.

## References

- Kuepper et al. (2015), posterior/MCMC and nuisance marginalization:
  https://arxiv.org/html/1502.02658
- Chen, Li & Gnedin (2025), measured progenitor MF and missing-mass correction:
  https://arxiv.org/pdf/2411.19899
- Chen, Gnedin & Price-Whelan, StarStream v2 (2026), Sections 4.1 and Appendix B:
  https://arxiv.org/html/2510.14924
- Baumgardt et al. (2023), MF convention and global MF measurements:
  https://arxiv.org/abs/2303.01636
- Official catalogue and estimated-slope legend:
  https://people.smp.uq.edu.au/HolgerBaumgardt/globular/parameter.html
