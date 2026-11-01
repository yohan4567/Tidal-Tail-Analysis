# Analysis allocation — 2026-11-01

Do not schedule new full numerical grids for the other-PC targets on the local
computer. Input preparation and downloads may continue locally. Scientific
holds are independent of which computer owns a target.

| Allocation | Targets |
|---|---|
| Local current run | M107 / NGC6171 |
| Local priority queue, in order | NGC6144, NGC5466, M5 / NGC5904, M80 / NGC6093, NGC2298, NGC3201 |
| Completed conditional analyses / preserved results | NGC288, Pal5 |
| Completed on other PC, user-reported; not imported locally | NGC362 |
| Other PC: ten-target conditional-input preparation | M4 / NGC6121, NGC6139, NGC6366, NGC6397, NGC6535, M62 / NGC6266, NGC6304, NGC6355, NGC6541, M2 / NGC7089 |
| Other PC: observation and model-input review | ESO452-SC11, NGC1851, M92 / NGC6341, M30 / NGC7099, Omega Cen / NGC5139, NGC5286, M19 / NGC6273, M22 / NGC6656, M15 / NGC7078 |

The allocation covers all 29 targets: 7 local, 3 completed/preserved and 19
assigned to the other PC. Completing a conditional analysis does not certify
an absolute physical total mass-loss rate.

| Review target | Current project-specific problem | Required next step |
|---|---|---|
| ESO452-SC11 | Saved global MF is indirect and has no measured mass interval; catalogue abundance/reddening fields need reconciliation. Age and composition literature exist. | Reconcile literature inputs; obtain a supported MF constraint or explicitly report population-dependent conditional results. |
| NGC1851 | Single-population conversion is unvalidated for its complex populations. | Review composition-dependent isochrones, population weights and their relevance to the escaping population. |
| M92 | Returned isochrone MH=-2.19174 differs from requested -2.31. | Obtain a supported composition-consistent model and validate the actual returned table. |
| M30 | Returned isochrone MH=-2.19174 differs from requested -2.27. | Same model check as M92; do not relabel the existing table. |
| Omega Cen | Multiple ages and compositions cannot be replaced by the catalogue's mean age. | Construct and check a weighted population mixture. |
| NGC5286 | Iron/s-process populations require review before an absolute star-count-to-mass conversion. | Review population-specific models and weights. |
| M19 | Known metallicity distribution is not represented by one isochrone. | Review population-specific models and weights. |
| M22 | Complex populations; abundance results depend on the spectroscopic analysis. | Review literature and CMD agreement before selecting a population mixture. |
| M15 | Existing MH=-2.19174 table lacks download metadata; the supplied file is not a verified replacement model. | Obtain a new supported model with provenance and parameter checks. |

References supporting the population review:

- [ESO452-11 age and composition](https://arxiv.org/abs/2502.01741).
- [NGC1851 multiple populations](https://doi.org/10.1051/0004-6361/202142234).
- [Omega Cen age/composition structure](https://arxiv.org/abs/2409.13855).
- [NGC5286 iron and s-elements](https://arxiv.org/abs/1502.07438).
- [M19 abundance distribution](https://arxiv.org/abs/1507.00756).
- [M22 differential abundance analysis](https://academic.oup.com/mnras/article/516/3/3515/6665944).

These review flags do not prohibit Gaia tail searches or orbit calculations.
They restrict the physical interpretation of an absolute mass-loss rate from
the current population conversion. The ten-target group also retains explicit
single-isochrone, escape-MF, extinction and calibration-transfer assumptions;
it is not automatically more scientifically validated than the review group.
