# Sparse quality calibration

NGC 6144 finished its 48 dynamics runs but failed during observation modelling: 9 of 40 quality cells had fewer than the required 20 eligible parent stars. All nine were in the blue G−RP bin. The original parent was a 10% independent Gaia sample; 152,823 stars remained after the recorded spatial exclusions.

The optional repair retains the original availability grid, magnitude bins, colour bins, and minimum eligible count. For each magnitude/colour stratum, if any spatial cell is sparse, all spatial cells in that stratum are pooled into a single quality-probability group. Other strata retain their original spatial resolution. Groups form a disjoint partition, so calibration counts are never duplicated across independent beta posteriors. Each group has one shared Beta(k+1,n−k+1) draw per calibration state.

This introduces an explicit assumption: the quality probability is spatially exchangeable within a pooled stratum. It supplies no new stars and does not validate the field-to-escape-population transfer. The grouping rule depends on eligible counts, not quality-success flags or a fitted escape rate. The current implementation does not remove unsupported stars or change the science fit mask.

The held-out field check subtracts held-out counts and successes before fitting. Predictions are checked in the original fine spatial cells, including the pooled strata. A maximum absolute beta-binomial predictive z of 3 is a declared diagnostic screen, not a calibrated goodness-of-fit test or evidence that the model is physically correct. Insufficient pooled or training support still raises an error.

## Separate, frozen experiment

```powershell
.\.venv\Scripts\python.exe repair_quality_calibration.py prepare --root outputs/cluster_survey_29_v1 --cluster NGC6144 --out outputs/NGC6144_quality_repair_v2
.\.venv\Scripts\python.exe -u outputs/NGC6144_quality_repair_v2/frozen_code/repair_quality_calibration.py validate --repair outputs/NGC6144_quality_repair_v2
```

Use another fresh repair folder for another experiment. Preparation records the selector location and source hashes and creates an isolated code snapshot. Only that snapshot loads the pooled map. The main production modules are unchanged, so a running queue continues with its original implementation. The representative validation retains 32 observation draws and 8 shared calibration states. It is one real-input response, not a complete grid or scientific recovery test.

After validation, run the snapshot:

```powershell
.\.venv\Scripts\python.exe -u outputs/NGC6144_quality_repair_v2/frozen_code/repair_quality_calibration.py run --repair outputs/NGC6144_quality_repair_v2 --workers 12
```

Existing dynamics are copied into the new grid only after numerical signatures, configuration, dependencies and file hashes pass the existing reuse checks. The old calibration, failed response attempt and dynamics remain unchanged. The three escape-MF offset assumptions remain separate. Results and status are stored in the repair folder, rather than replacing the original survey entry.

To avoid competing with an existing queue, the Windows wrapper can wait for that queue's actual PID:

```powershell
powershell -NoProfile -File scripts/run_quality_repair_after_queue.ps1 -RepairPath outputs/NGC6144_quality_repair_v2 -QueuePid 54768 -Workers 12
```

Replace the PID with the verified active queue PID. `scheduled_run.json` reports waiting/running/failed/completed states; `full_run.log` records the numerical run. This wrapper must itself remain running. It does not automatically relaunch after shutdown, and only one writer may use a repair folder. If interrupted, inspect its status and rerun the frozen `run` command when no other writer remains.

## Local checks on 1 November 2026

For the recorded NGC 6144 inputs, five blue magnitude/colour strata were pooled, producing 25 quality groups instead of 40 independent cells. The smallest group had 54 eligible stars. The original-fine-cell holdout diagnostic tested 38 cells and had maximum absolute predictive z = 1.8795, with none above 3. The representative response completed without the original low-count-cell error. Full-grid recalculation and scientific validation are separate subsequent stages; no validated physical escape rate follows from these checks alone.
