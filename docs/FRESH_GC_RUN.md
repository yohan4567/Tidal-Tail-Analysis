# Fresh analysis: NGC 6093 (M80), NGC 1851, and NGC 6341 (M92)

Run these commands in PowerShell on the destination PC, from the repository directory. This starts new observations and simulations without copying old checkpoints. Keep the new output directory local to that PC, with only one queue writing to it.

## Requirements and working directory

Install Git and 64-bit Python 3.11 on the destination PC. Check them in PowerShell:

```powershell
git --version
py -3.11 --version
```

The examples use `C:\GCResearch`, outside OneDrive. Create it and start there:

```powershell
New-Item -ItemType Directory -Force C:\GCResearch | Out-Null
Set-Location C:\GCResearch
```

For the 12-worker, 24-GiB soft-budget example below, use a machine with sufficient CPU capacity and more than 24 GiB of installed RAM; 32 GiB is the development configuration. Other active applications also consume RAM. Download sizes and runtimes depend on the real field and calibration downloads; no completion time is guaranteed.

## Install

```powershell
git clone https://github.com/yohan4567/Tidal-Tail-Analysis.git
cd Tidal-Tail-Analysis
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
if ($LASTEXITCODE -ne 0) { throw 'Dependency installation failed.' }
.\.venv\Scripts\python.exe -m pip install 'requests>=2.31,<3'
if ($LASTEXITCODE -ne 0) { throw 'Downloader dependency installation failed.' }
.\.venv\Scripts\python.exe -m pip check
```

For an existing checkout, update it with `git pull --ff-only` before starting instead of cloning again.

`requests` is required by the external field downloader; it is not a core dependency of this repository. All Python commands below use the virtual environment directly, so PowerShell activation and execution-policy changes are unnecessary.

## Prepare real inputs

The external `tidal_study_v2` project and its data are not included in this repository. Place an authorized copy on the destination PC and set its path:

```powershell
$env:GC_TIDAL_STUDY = 'C:\GCResearch\tidal_study_v2'
```

Replace the example path with the actual directory containing `step3_tail_search.py`. It must also contain its imported modules, `catalogs/master_clusters.csv`, and each cluster's real Gaia field and matching query metadata. Missing-tile masks, where applicable, must accompany the fields. Copying code alone does not supply these observations. Set this environment variable again in each new PowerShell session.

The commands below deliberately use the original field layout (`--legacy`). Use it only if those are the fields you intend to analyze. To use prepared v3 fields instead, remove `--legacy`; the software does not silently switch between them.

The expected original-field layout is:

```text
C:\GCResearch\
  Tidal-Tail-Analysis\
    .venv\
    run_pipeline.py
  tidal_study_v2\
    step2_query_gaia_fields.py
    step3_tail_search.py
    skymath.py
    catalogs\master_clusters.csv
    gaia_fields\
      NGC6093_field.csv.gz
      NGC6093_field.json
      NGC1851_field.csv.gz
      NGC1851_field.json
      NGC6341_field.csv.gz
      NGC6341_field.json
```

Keep the rest of the external project too. The tree only highlights the relevant files. Transfer corresponding `*_missing_tiles.json` files whenever present. Transfer raw inputs, not the old virtual environment or old simulation outputs.

### Download missing observed fields

Skip this subsection if all three original fields and metadata are already available. The following commands require the inspected version of the external downloader; check its help/documentation if using a different version. Remain in `C:\GCResearch\Tidal-Tail-Analysis` so `.venv` resolves correctly.

The normal target list includes all three selected clusters:

```powershell
.\.venv\Scripts\python.exe "$env:GC_TIDAL_STUDY\step2_query_gaia_fields.py" --release=dr3 NGC6093 NGC1851 NGC6341
```

Failed downloads are logged even when the downloader exits zero. Check the files below, not just the exit code. Keep matching CSV files and query metadata together.

On repeated request timeouts, retry the same failed target with smaller chunks (`--tiles=2`); `--level=8` additionally changes the sky tiling and query-plan metadata. Treat a changed query plan as a new input preparation, before freezing any analysis. Do not mix data, query metadata, and missing-tile masks from different download plans. Avoid `--skip-failed` for an initial complete-field attempt.

### Check inputs before analysis

```powershell
@'
from pathlib import Path
import json
import pandas as pd
from selection_bridge import DEFAULT_STUDY, REQUIRED, load_study

master = pd.read_csv(DEFAULT_STUDY / 'catalogs/master_clusters.csv').set_index('id')
selector = load_study(DEFAULT_STUDY, legacy=True)
for cid in ('NGC6093', 'NGC1851', 'NGC6341'):
    if cid not in master.index:
        raise ValueError(f'Missing master-catalogue row: {cid}')
    field = selector.field_file(cid, '_field.csv.gz')
    metadata = selector.field_file(cid, '_field.json')
    if not field.is_file() or not metadata.is_file():
        raise FileNotFoundError(f'{cid}: need both {field} and {metadata}')
    plan = json.loads(metadata.read_text(encoding='utf-8'))
    if not all(k in plan for k in ('radius_deg', 'pm0', 'pm_box', 'g_limit')):
        raise ValueError(f'Incomplete query metadata: {metadata}')
    sample = pd.read_csv(field, nrows=5)
    missing = set(REQUIRED) - set(sample.columns)
    if missing:
        raise ValueError(f'{cid}: missing Gaia columns {sorted(missing)}')
    if sample.empty:
        raise ValueError(f'Empty Gaia field: {cid}')
    print(f'{cid}: input files and columns found; radius={plan["radius_deg"]} deg; {field}')
'@ | .\.venv\Scripts\python.exe -
if ($LASTEXITCODE -ne 0) { throw 'Fix the input errors before continuing.' }
```

This checks file presence, metadata keys, and sample columns. The full pipeline performs additional consistency checks; this is not a completeness or scientific validation test.

## Analyze the three observed fields

Choose a new folder name for a genuinely new experiment. Existing pipeline output directories are not overwritten.

```powershell
$root = 'outputs/new_gcs_6093_1851_6341'
$clusters = @('NGC6093', 'NGC1851', 'NGC6341')

foreach ($cluster in $clusters) {
    .\.venv\Scripts\python.exe run_pipeline.py --cluster $cluster --legacy --out "$root/$cluster"
    if ($LASTEXITCODE -ne 0) {
        Write-Warning "Observed analysis failed for $cluster; check its logs. Continuing."
        continue
    }
    .\.venv\Scripts\python.exe make_config.py $cluster --out "$root/$cluster/adopted_dynamics.json"
    if ($LASTEXITCODE -ne 0) {
        Write-Warning "Dynamics configuration failed for $cluster; inspect the error before inference."
    }
}
```

Each successful cluster has a `REPORT.md`, `selection/selection.json`, selection maps, and a `tail_extent/` directory. These are conditional observed tail-candidate profiles; this stage does not infer a mass-loss rate.

## Prepare shared inference data

Download the real completeness map and the measured mass-function source table once. Existing files are retained; source hashes are recorded during subsequent preparation.

```powershell
if (-not (Test-Path 'outputs/cache/gaia/allsky_uniq_10.fits')) {
    .\.venv\Scripts\python.exe selection_model.py fetch-gaia-map --mode multi --out outputs/cache/gaia/allsky_uniq_10.fits
    if ($LASTEXITCODE -ne 0) { throw 'Completeness-map download failed.' }
}
if (-not (Test-Path 'outputs/empirical_mf_source/parameter.html')) {
    New-Item -ItemType Directory -Force outputs/empirical_mf_source | Out-Null
    Invoke-WebRequest 'https://people.smp.uq.edu.au/HolgerBaumgardt/globular/parameter.html' -OutFile outputs/empirical_mf_source/parameter.html -ErrorAction Stop
}
```

## Run conditional mass-loss inference

```powershell
.\.venv\Scripts\python.exe -u run_deadline_queue.py --root $root --clusters NGC6093 NGC1851 NGC6341 --workers 12 --memory-budget-gib 24
```

Use 12 workers only on a suitable PC; reduce `--workers` and the RAM budget for a smaller machine. The RAM budget is a soft scheduling/headroom check, not an operating-system memory limit. The queue runs one cluster at a time and uses workers within its numerical stages. It prepares missing official isochrones and Gaia calibration samples, so network access is needed. It logs a failed or blocked cluster and continues to the next.

NGC 6093 can enter conditional inference once its real inputs and calibrations pass preparation. Its supplemental population entry records the scaled-solar abundance approximation used by the baseline workflow.

NGC 6341 is configured as a target, but an earlier local preparation failed the isochrone-metallicity consistency check: the adopted catalogue value is [Fe/H] = -2.31, whereas the returned model did not match the declared metallicity. A fresh download is checked again; the queue may stop preparation for this cluster until a supported, composition-consistent population model is supplied. Do not relabel a clipped isochrone or bypass that check to force a result. Its observed analysis is independent of this inference blocker.

Successfully finishing a conditional run does not establish a validated absolute physical mass-loss rate.

NGC 1851 retains a population-review hold. Its observed analysis can run, but the mass-loss queue records it as blocked because the current single-population model has not been approved for that cluster. Do not remove the hold merely to obtain a number.

## Monitor and resume

In a second PowerShell window, from the repository directory:

```powershell
Get-Content outputs/new_gcs_6093_1851_6341/EXECUTION_LOG.md
Get-Content outputs/new_gcs_6093_1851_6341/execution_events.jsonl -Tail 10 -Wait
```

The Markdown file summarizes the queue stages. The JSONL stream records commands, stage changes, and resource samples; it is not a particle-level percentage display. Press Ctrl+C to stop watching the log.

After an orderly interruption, rerun the same queue command to reuse valid checkpoints. A recorded failure is skipped unless `--retry-failed` is supplied after its cause is corrected. Do not rerun the observed pipeline into an existing cluster directory; retain its successful selection when resuming inference. Results for successful inference live under `<root>/<cluster>/empirical_direct_v1/`, with separate offset scenarios.

For example, after fixing a failed download or calibration input:

```powershell
.\.venv\Scripts\python.exe -u run_deadline_queue.py --root outputs/new_gcs_6093_1851_6341 --clusters NGC6093 NGC1851 NGC6341 --workers 12 --memory-budget-gib 24 --retry-failed
```

An expected NGC 1851 population hold is not fixed by retrying. For a new PowerShell session, first return to the repository and restore the external-project environment variables. The literal output path in the resume command avoids relying on the old `$root` variable.

Keep the PC connected to power and prevent automatic sleep during computation. Closing the terminal or shutting down can interrupt the queue; after shutdown, manually restart it. Checkpoints preserve completed work where supported, but an interrupted in-flight chunk may need recomputation. Never start a second writer in the same output directory.

## Read the results

The top-level observed `REPORT.md` is not updated into a mass-loss report after queue completion. Look under the inference directory for those products:

```text
outputs/new_gcs_6093_1851_6341/
  EXECUTION_LOG.md
  execution_events.jsonl
  deadline_queue_status.json
  NGC6093/
    REPORT.md
    selection/
    tail_extent/
    adopted_dynamics.json
    calibration_inputs/
    science_calibration/
    empirical_direct_v1/
      transfer_sensitivity.json
      offset_0/posterior/
      offset_minus0p5/posterior/
      offset_minus1/posterior/
```

Each successful posterior directory contains `posterior_summary.json`, `joint_posterior.png`, and `posterior_prior_sensitivity.png`. The three offsets are separate escape-population assumptions; do not average them as if they were measured equal-probability alternatives. A queue status of `finished` means every requested target was processed, including blocked or failed targets, not that every target produced a valid posterior. Inspect individual entries in `deadline_queue_status.json`.

Common issues:

| Message or symptom | Action |
|---|---|
| External selector missing | Point `GC_TIDAL_STUDY` at the directory containing `step3_tail_search.py`. |
| Missing observed field or metadata | Restore/download the matching pair and any missing-tile mask; check the resolved paths printed by the input check. |
| Output directory already exists | Use another root for a new observation analysis; use the existing root only to resume its inference. |
| Missing `adopted_dynamics.json` | Run `make_config.py` for that cluster after its observed stage succeeds. |
| TCG map or quality grid has no support | Investigate data/calibration coverage. Do not replace unsupported probabilities with invented values. |
| Frozen-input checksum mismatch | Restore the frozen input, or start a separate experiment after documenting the change. |
| Insufficient RAM headroom | Close other memory-heavy programs; choose workers and a soft budget suitable for installed RAM. Increasing the budget does not add physical RAM. |
| galpy C-extension warning | The project has a tested Cartesian CPU fallback; inspect the actual final error/status before treating a warning as a failure. |

The present conditional baseline leaves population transfer, abundance mapping, background and potential assumptions to be assessed. Passing setup checks or obtaining a posterior is distinct from completing the scientific validation needed to claim an absolute physical mass-loss rate.
