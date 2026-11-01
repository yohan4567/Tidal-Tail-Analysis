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
use their new versioned jobs and retain the old failed outputs. A completed retry
is excluded from the later general queue to avoid rerunning the broken inputs.

### Compact-path recovery

Windows cache paths in the first preparation exceeded the normal path-length
limit. `recover_compact_repairs.py` creates separate compact preparations and
uses `outputs/rc2/<cluster>/<hash>` for caches. Path-only repairs are compared
against the original real-response arrays, including all calibration draws.
Physical observation settings, sample counts, and donor distances are unchanged.

M80 additionally uses the existing spatial pooling within each magnitude/colour
stratum when its original spatial cells have insufficient calibration support.
The minimum count stays 20. Pooling requires a held-out spatial check and retains
a common Beta draw for cells in the same group. Availability calibration retains
its original spatial grid. This assumes exchangeability within the pooled
stratum; passing the check is not proof of full spatial calibration. A response
using the existing full-duration simulation is required before the local retry.

`prepare_exact_local_retries.py` creates new M5/M80 jobs with the original
dynamics configuration. It permits an annotation-only difference in
`mass_evolution.source_description`; any physical or numerical difference is
rejected. The strict grid-reuse guard is not relaxed. Old failed runs are kept.

`resume_indexed_calibrations.py` retries M62 and other failed download targets
with level-7 HEALPix patches, one per request. Completed matching patches remain
reusable. No partial download or ZIP-only export qualifies as ready.
`prepare_ngc3201_bright_retry.py` separately restores the missing bright-census
support for NGC 3201 and waits behind M5/M80 before its optional local run.

### Appended local queue

The subsequent user instruction assigns the seven repaired targets to the local
queue. `queue_repaired_seven.py` waits for the existing M5/M80 dispatcher and
NGC 3201 retry to finish. It records M5/M80's already completed results rather
than rerunning them, then processes M92, M30, M15, NGC 6139, and M62 serially.
M62 waits for its validated package; a recorded preparation or analysis failure
is retained and the queue continues. The other targets' portable packages stay
in place for transfer. Previously reserved targets outside these seven retain
their allocation.

For the November 1 recovery roots:

```powershell
.\.venv\Scripts\python.exe -u queue_repaired_seven.py --root outputs/q7r --workers 12
```

Only one instance may own an output root. `plan.json` records the immutable queue
configuration; `status.json`, `queue.log`, and per-cluster logs record progress.
The queue checks archive and manifest identities against the relocation test
before starting each job. Completed and failed entries remain terminal on a
restart; changing a failed model requires a reviewed new job/output root. Twelve
workers belong to one cluster at a time. This command does not start another
grid while the specified predecessor jobs are active.

Progress is recorded in `status.json`, `prepare.log`, and the per-cluster logs.
Only packages with a passed relocated real-input response and a verified ZIP CRC
are added to this run's `READY_PACKAGES.json`. Export alone is not readiness.
New exports use likelihood comparisons, with no R/T posterior. The actual
isochrone age and explicit population assumptions travel with the package.

Passing these software and input checks permits a conditional analysis. Full
template convergence, recovery/coverage tests, population/selection calibration,
and physical systematic checks remain necessary for an absolute-rate claim.
