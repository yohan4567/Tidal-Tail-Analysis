# Versioned survey input repairs

These repairs preserve the existing survey outputs and numerical integration
sources. Run them in a new output directory; do not overwrite a running grid.
The repair set is NGC 3201, NGC 362, M4 (NGC 6121), M107 (NGC 6171), NGC 6397,
M62 (NGC 6266), and M2 (NGC 7089).

## What changes

- Rounded main-sequence current-mass maps can assign several photometric nodes
  the same mass coordinate. Integrate the power-law measure on unique mass
  coordinates, then allocate each tied group's measure using the original
  isochrone phase weights. This retains every photometric node, conserves the
  group's measure, and normalizes to one unit of represented current stellar
  mass. The phase allocation is an explicit approximation, recorded in metadata;
  it does not measure the escaping mass function or extend the isochrone.
- Parent downloads split a batch of sky tiles on TAP `OVERFLOW`. A single tile
  is refined into its four nested HEALPix children if necessary. Only complete
  responses are accepted. Cached queries must match exactly and pass checksum
  checks; duplicate or out-of-tile source IDs cause failure. Sky, G and independent
  random-index sampling are retained, with no membership or quality cuts.
- Checkpoint replacement retries `PermissionError` up to 12 times, with bounded
  backoff. Persistent errors still stop that cluster. The I/O wrapper leaves the
  integration sources, configurations and checkpoint identities unchanged.
  Existing simulations/checkpoints are reused only after exact source, config,
  dependency and checksum checks.

The root production modules are left unchanged while the original queue is
active. `repair_survey_inputs.py bootstrap` creates a frozen source snapshot with
the changes above and hashes its inputs. The repair run uses this snapshot.

## Commands (PowerShell, from the repository root)

The original survey selection, adoption files, official isochrones, global MF
table and official Gaia completeness map must already be available locally.
Set `GC_TIDAL_STUDY` to the matching external selection repository if needed.

```powershell
$repair = 'outputs/input_repairs_new_version'
.\.venv\Scripts\python.exe repair_survey_inputs.py bootstrap --source outputs/cluster_survey_29_v1 --out $repair
.\.venv\Scripts\python.exe "$repair/frozen_code/repair_survey_inputs.py" prepare --out $repair
# Run only after other full analysis queues have stopped:
.\.venv\Scripts\python.exe "$repair/frozen_code/repair_survey_inputs.py" run --out $repair --workers 12
```

Input preparation is serial and continues after a cluster fails. Its status is
`input_status.json`; downloads also update each cluster's
`calibration_inputs/adaptive_parent_status.json`. Read `prepare.log` when launched
in the background. After preparation, the run command selects only ready input
sets, continues after cluster failures, and writes per-cluster
`inference_status.json`. `ready_for_conditional_inference` means input construction
checks passed; response support and independent physical validation remain.

`scripts/run_input_repairs_after_queue.ps1` can wait for specified existing job
PIDs and the preparation PID before starting the repaired clusters, avoiding two
12-worker analysis queues at once. It records `scheduled_run.json`. A finished
scheduler does not imply every cluster succeeded; inspect the individual results.

These repairs do not certify an absolute mass-loss rate. Selection transfer,
escape-population assumptions, dynamical applicability, numerical/Monte Carlo
convergence, and recovery/coverage tests remain separate requirements.

## Verification

```powershell
.\.venv\Scripts\python.exe -m unittest test_input_repairs test_batch_dynamics -v
```

Tests check analytic mass-function integrals, tied-node measure and mass
normalization, exact recursive sky coverage/cache integrity, bounded failures,
and real integration/checkpoint resumption after an injected file lock.
On managed Windows systems, run file-replacement tests with normal filesystem
access: sandbox denial of replacement is persistent and cannot be fixed by retrying.
