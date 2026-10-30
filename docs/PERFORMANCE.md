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
workers, batch size 128 and two response workers. It sets BLAS/OpenMP thread
limits to one per process to avoid nested oversubscription. These are worker
processes, not hard affinity assignments to particular physical CPU cores.
Data/calibration inputs are external and are not included in a clean clone.

```powershell
.\scripts\Run-Pal5Fast.ps1 -PlanOnly
.\scripts\Run-Pal5Fast.ps1
```

The live accelerated output is
`outputs/cluster_survey_29_v1/Pal5/production_grid_fast_v1`. Only rerun after the
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
