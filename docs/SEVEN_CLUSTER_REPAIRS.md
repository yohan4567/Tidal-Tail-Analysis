# Seven-cluster input repair

The requested targets are M92 (NGC 6341), M30 (NGC 7099), M15 (NGC 7078),
M5 (NGC 5904), M62 (NGC 6266), NGC 6139, and M80 (NGC 6093).
All repairs use separate output directories. Existing completed analyses,
source-frozen snapshots, selection masks, and packages remain intact.

## Composition inputs

`configs/four_population_repairs_v1.json` records the adopted literature values
and their limitations. For M92, M30, and M15, the age study's iron abundance and
adopted alpha enhancement are converted into an explicitly approximate global
metallicity for scaled-solar PARSEC tracks. This is not a directly alpha-enhanced
stellar model. The returned metallicity and age must match the request; clipped
models are rejected. A separate alpha=0.30 sensitivity isochrone is retained,
not treated as a measured confidence bound. The nominal model drives the run.

The original selection and distance remain fixed. CMD ridge comparisons and
positive, unit-current-stellar-mass population quadratures are saved for review.
These checks do not validate horizontal-branch weights, differential extinction,
the escaping mass function, remnants, or an absolute physical escape rate.

## Gaia calibration repair

`repair_parent_support.py` tries bounded, complete Gaia TAP queries, saves endpoint
errors, and switches endpoints on transport/HTTP failures. VOTable overflow is
never accepted as a complete parent. Larger faint samples use disjoint adaptive
HEALPix queries; completed components are verified by checksum before reuse.
The merge can wait for RAM without losing downloaded tiles.
An HTTP 400 whose VOTable explicitly reports a statement timeout triggers indexed
tile queries instead of treating the query as a successful or merely truncated
result. A timed-out indexed patch is subdivided further.

Calibration uses an uncut 10% random-index sample at 12 <= G < 21 and a separate
complete bright census at 0 <= G < 12. Components must have the same sky geometry,
documented independent sampling, and no membership or quality cuts. The existing
science field is not substituted for a calibration parent.

M80's failure is checked again using an existing full-duration simulation after
the calibration is rebuilt. A 16-particle smoke test alone cannot establish that
the original full-duration donor-support problem has been repaired. Support
thresholds and donor feature distances are not silently relaxed.

## Execution and readiness

From the repository root, first prepare the three composition inputs:

```powershell
.\.venv\Scripts\python.exe -u repair_four_populations.py --root outputs/seven_cluster_repairs_20261101_v1/composition --clusters NGC6341 NGC7099 NGC7078 --inputs-only
```

After that command completes, prepare and validate the seven targets:

```powershell
.\.venv\Scripts\python.exe -u repair_seven_clusters.py --root outputs/seven_cluster_repairs_20261101_v1
```

Use one writer per output directory. The worker downloads inputs before running
memory-intensive calibrations, waits for at least 8 GiB of available RAM, and
continues to the next target after a recorded failure. Existing numerical
analyses are not stopped. This preparation command does not launch concurrent
12-worker production grids.

An explicit `ready_local_repairs_root` setting in the old repair dispatch
overrides lets the existing scheduler retry M5 and M80 after its current priority
pass, serially and only after the repaired packages pass validation. These retries
use their new versioned jobs and retain the old failed outputs. Other targets
remain reserved for the other computer. Only a completed retry is excluded from
the later general queue to avoid rerunning the original broken configuration.

Progress is recorded in `status.json`, `prepare.log`, and the per-cluster logs.
Only packages with a passed relocated real-input response and a verified ZIP CRC
are added to this run's `READY_PACKAGES.json`. Export alone is not readiness.
New exports use likelihood comparisons, with no R/T posterior. The actual
isochrone age and explicit population assumptions travel with the package.

Passing these software and input checks permits a conditional analysis. Full
template convergence, recovery/coverage tests, population/selection calibration,
and physical systematic checks remain necessary for an absolute-rate claim.
