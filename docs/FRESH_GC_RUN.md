# Fresh analysis: NGC 6139, NGC 6093 (M80), and NGC 1851

Run these commands in PowerShell on the destination PC, from the repository directory. This starts new observations and simulations without copying old checkpoints. Keep the new output directory local to that PC, with only one queue writing to it.

## Install

```powershell
git clone https://github.com/yohan4567/Tidal-Tail-Analysis.git
cd Tidal-Tail-Analysis
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

For an existing checkout, update it with `git pull --ff-only` before starting instead of cloning again.

## Prepare real inputs

The external `tidal_study_v2` project and its data are not included in this repository. Place an authorized copy on the destination PC and set its path:

```powershell
$env:GC_TIDAL_STUDY = 'D:\research\tidal_study_v2'
```

Replace the example path with the actual directory containing `step3_tail_search.py`. It must also contain its imported modules, `catalogs/master_clusters.csv`, and each cluster's real Gaia field and matching query metadata. Missing-tile masks, where applicable, must accompany the fields. Copying code alone does not supply these observations. Set this environment variable again in each new PowerShell session.

The commands below deliberately use the original field layout (`--legacy`). Use it only if those are the fields you intend to analyze. To use prepared v3 fields instead, remove `--legacy`; the software does not silently switch between them.

## Analyze the three observed fields

Choose a new folder name for a genuinely new experiment. Existing pipeline output directories are not overwritten.

```powershell
$root = 'outputs/new_gcs_6139_6093_1851'
$clusters = @('NGC6139', 'NGC6093', 'NGC1851')

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
.\.venv\Scripts\python.exe -u run_deadline_queue.py --root $root --clusters NGC6139 NGC6093 NGC1851 --workers 12 --memory-budget-gib 24
```

Use 12 workers only on a suitable PC; reduce `--workers` and the RAM budget for a smaller machine. The RAM budget is a soft scheduling/headroom check, not an operating-system memory limit. The queue runs one cluster at a time and uses workers within its numerical stages. It prepares missing official isochrones and Gaia calibration samples, so network access is needed. It logs a failed or blocked cluster and continues to the next.

The current configuration allows conditional inference for NGC 6139 and NGC 6093. Their supplemental population entries explicitly record the same scaled-solar abundance approximation used by the baseline workflow. NGC 6139's differential reddening remains unmodelled. Successfully finishing these runs does not establish a validated absolute physical mass-loss rate.

NGC 1851 retains a population-review hold. Its observed analysis can run, but the mass-loss queue records it as blocked because the current single-population model has not been approved for that cluster. Do not remove the hold merely to obtain a number.

## Monitor and resume

In a second PowerShell window, from the repository directory:

```powershell
Get-Content outputs/new_gcs_6139_6093_1851/EXECUTION_LOG.md
Get-Content outputs/new_gcs_6139_6093_1851/execution_events.jsonl -Tail 10 -Wait
```

The Markdown file summarizes the queue stages. The JSONL stream records commands, stage changes, and resource samples; it is not a particle-level percentage display. Press Ctrl+C to stop watching the log.

After an orderly interruption, rerun the same queue command to reuse valid checkpoints. A recorded failure is skipped unless `--retry-failed` is supplied after its cause is corrected. Do not rerun the observed pipeline into an existing cluster directory; retain its successful selection when resuming inference. Results for successful inference live under `<root>/<cluster>/empirical_direct_v1/`, with separate offset scenarios.
