# Batched integration and CPU parallelism

The `cartesian_batch` backend evaluates independent particle forces together in
NumPy while using the same MWPotential2014 components, moving Plummer potential
and release-state sampler. It uses the incomplete-gamma expression equivalent
to galpy's cutoff-bulge enclosed mass, and combines shared spherical force terms.
A 1000-position test spanning Galactic radii compares forces to galpy at 1e-10
relative tolerance. Unsupported potential component classes are rejected.

Particles with different release times use a normalized integration coordinate
u in [0,1], physical time t=age*(u-1). The RHS includes the age Jacobian. A DOP853
subclass takes the maximum error norm over individual six-coordinate tracers;
adding well-behaved tracers therefore cannot dilute one poorly integrated tracer.
The protected SciPy callback is covered by regression tests; SciPy version is
recorded in the simulation identity and metadata.

Production uses rtol=1e-11 and atol=1e-13, tighter than the original scalar
1e-9/1e-11 settings. Backend and batch size are explicit configuration fields.
Changing the backend changes checkpoint identity. A completed batch is saved to
an atomic checkpoint once at least 100 new tracers have accumulated. For the
128-particle setting this saves every batch: abrupt interruption loses at most
the unfinished batch, rather than the entire simulation. Configuration, source
hashes and package versions must match to resume; never start two writers in the
same output directory.

## Measured checks on Pal 5

A representative 128-tracer sample from nearby release-time strata at about
1.5 Gyr was compared with a stricter scalar reference (rtol=1e-12,
atol=1e-14). The batch-128 run took 2.60 seconds versus 219.99 seconds for
that stricter reference. Maximum differences were 1.53e-5 pc and 1.86e-7 km/s.
The original scalar tolerance took 63.83 seconds on the same layout in an earlier
run; comparison with 2.60 seconds suggests about 24.6x improvement on this small
sample, not a guaranteed whole-pipeline speedup. Record CPU time as well as wall
time because sleep, contention and scheduling affect wall time.

The first adjacent-stratum comparison at the original tolerance exceeded the
prespecified 0.01 pc difference gate (about 0.013 pc). It is retained as a failed
numerical comparison; the gate was not relaxed. Production uses the tighter
settings that passed the stricter-reference comparison. End-to-end 96-tracer
simulations at rates 2, 9 and 32 were each audited on eight release states
spanning the modeled window; all passed 0.01 pc and 0.001 km/s gates. These checks
are numerical checks, not validation of population assumptions or physical MLR.

## Ten-worker run and restart

`scripts/Run-Pal5Fast.ps1` runs the prepared local Pal 5 inputs with 10 dynamics
workers, batch size 128 and six response workers. Two response workers were
initially used. Ten were tested after measuring about 280 MB resident memory per
process, but whole-system available RAM subsequently fell below 1 GB. Six were
selected to leave more headroom. Resident memory alone does not capture peak or
whole-system demand. Check available RAM and throughput before increasing this
setting on another machine or workload. It sets BLAS/OpenMP thread
limits to one per process to avoid nested oversubscription. These are worker
processes, not hard affinity assignments to particular physical CPU cores.
Data/calibration inputs are external and are not included in a clean clone.

```powershell
.\scripts\Run-Pal5Fast.ps1 -PlanOnly
.\scripts\Run-Pal5Fast.ps1
```

The live accelerated output is
`outputs/cluster_survey_29_v1/Pal5/production_grid_odd_v1`. It uses the requested
16 rates (1, 3, ..., 31) and 20 time windows: 320 nodes and 48 dynamics runs.
Matching odd-rate work from `production_grid_fast_v1` is reused after strict
provenance checks; old results remain intact. Only rerun after the
previous writer has stopped. Normal shutdown preserves saved files but does not
automatically restart the job when the computer boots. Rerunning with identical
inputs/code resumes completed jobs and integration checkpoints. Updates to code
require retaining the matching frozen runtime or starting a fresh output.

The old runtime is archived locally in `outputs/Pal5_reference_runtime_v1` and
old results remain in `outputs/nine_cluster_real_inputs_v1/Pal5/production_grid`.
Because the accelerated run uses tighter tolerance, reference checkpoints are
not silently merged into the new run. The survey progress pointer identifies
the new output; the visible progress count consequently restarts for the new
numerical experiment.

## Grid cost

Each extra rate adds three independent dynamics runs at the current seed count.
At Tmax=2000 Myr and a 0.25 Myr per-arm interval, that is 48000 extra trajectories.
Extra T windows within an unchanged Tmax reuse trajectories and mostly add
response-window bookkeeping and likelihood scores. Extending Tmax increases
both particle count and mean integration age; a roughly quadratic cost is a
useful heuristic at fixed release interval, not an exact timing law.

Exhaustive per-truth confidence calibration still scales roughly as K squared
in candidate likelihood fits. The batched dynamics backend does not speed up
that statistical stage. Data-dependent coarse-to-fine search must be included
in recovery experiments if a confidence claim depends on the search procedure.


The response concurrency change preserves all scientific code, random seeds,
population/calibration draws and completed templates. The two interrupted
in-flight response calculations are rerun if no complete output exists. No
sample-count reduction or approximate replacement is involved. Higher
concurrency does not guarantee proportional speedup; CPU, RAM and disk
contention must be measured after restart.


## Shared observations (production default)

`response_model.py` now shares observation draws across quality states by default.
`response_group.py` also shares those draws across compatible population variants,
and `run_joint_study.py` schedules one response group per dynamics realization.
All population nodes, calibration states and observation draws are retained. Each
population keeps its own native response files, MC standard errors, cutflow and
per-tracer contributions. Numerical backend names and source hashes are recorded
in response provenance, the cache key and the grid plan.

Sharing across populations requires identical ordered nodes, photometry, mass basis,
positive-weight support, observation settings and calibration inputs. Incompatible
populations/settings are processed separately with sharing within their quality
states. Different dynamics/MLR/seeds are never combined. The original reference
implementation remains available with `response_model.py --backend reference`.

Validation on Pal 5 real calibration inputs:

| Particle subset | Cases | Noise draws | Reference | Shared | Kernel speedup |
|---|---:|---:|---:|---:|---:|
| 128 | 24 | 32 | 12.925 s | 1.053 s | 12.27x |
| 1,024 | 24 | 32 | 96.978 s | 6.253 s | 15.51x |

These are single-machine subset timings, not end-to-end runtime claims. The full
16,000-particle production export was separately compared against previously
completed reference outputs for three populations: **every array** in calibrated
counts, response maps (including MC errors) and tracer responses was bitwise equal.
Age-contribution CSVs, cutflows, sampling identities, 100/1,000 Myr derived templates
and completed-group resume also matched or passed validation. These numerical
checks do not validate the underlying physical assumptions or absolute MLR.

`benchmark_response_shared.py` reproduces the kernel comparison on existing real
inputs and writes source hashes, equality checks and timings. For example:

```powershell
$grid = 'outputs/cluster_survey_29_v1/Pal5/production_grid_odd_v1'
$configs = 0..2 | ForEach-Object { "$grid/configs/rate_9_seed_1009_pop_$($_)_T_2000.json" }
.\.venv\Scripts\python.exe benchmark_response_shared.py --selection outputs/cluster_survey_29_v1/Pal5/selection --simulation "$grid/rate_9_seed_1009/simulation" --configs $configs --particles 1024 --output outputs/shared_response_benchmark_1024.json
```

The Pal 5 launcher uses `production_grid_shared_v1`, three response-group workers
and ten dynamics workers. Its inactive `production_grid_odd_v1` source retains the
original outputs. All 48 existing dynamics can be reused after exact checks, while
response templates are rebuilt under the new backend/provenance. An old frozen
plan is never rewritten to bypass a source mismatch. Monitor actual free memory
before raising group concurrency; one group now represents up to 24 cases.
Completed valid response groups are skipped on resume. Incomplete output folders
are refused for inspection rather than overwritten. Do not run two writers
against the same output grid.

## Direct empirical queue (2026-10-31)

`run_remaining_survey.py --prepare-only` prepares cluster-specific inputs without
executing a legacy tilt grid. `run_direct_empirical.py --cluster NGC5466 --workers 12`
then computes the measured-MF models directly: 16 rates, 20 time windows, three
dynamics seeds, 16,000 tracers per dynamics realization, five population nodes
and three separate escape-transfer offsets. These numerical settings are unchanged.
The first offset supplies verified dynamics reuse for subsequent offsets. Existing
inactive legacy dynamics can also be reused if exact compatibility checks pass.

This omits the extra three-population legacy response grid and its posterior
post-processing. It does not claim a measured end-to-end speedup or new physical
validation. Twelve independent subprocesses handle templates/dynamics; background
quadrature remains serial. BLAS/OpenMP threads are limited to one per subprocess.

`run_deadline_queue.py --workers 12 --memory-budget-gib 24` executes one cluster at
a time, prioritizing interrupted work and then prepared inputs. Existing completed
posteriors are skipped; Pal5's cancelled rerun remains preserved. Previously failed
clusters are skipped unless `--retry-failed` is explicitly supplied after correction.
Blocked input populations are recorded, never replaced by another cluster's data.
The queue continues after an individual failure. Progress is recorded in
`outputs/cluster_survey_29_v1/deadline_queue_status.json`, with per-cluster
`deadline_inference.log` and `empirical_job.json`.
`scripts/run_continuous_survey.ps1` provides the Windows launcher; it does not wait
for a historical PID. The 24 GiB setting is a soft system-headroom policy, not an
OS memory cap or an artificial allocation. Commands, exit codes, errors and regular
system-RAM samples are appended to `execution_events.jsonl`; `EXECUTION_LOG.md`
summarizes the current pass. No integration batches, grids or seeds change.
Do not start a second writer against the same directories.
Completion time depends on data availability; this queue cannot guarantee all
29 clusters finish by a deadline. Scientific validation is deferred, not passed.
