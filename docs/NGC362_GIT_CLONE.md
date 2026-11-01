# NGC 362: Git clone setup on another Windows PC

Use PowerShell. Git supplies the project checkout; the separately transferred
`NGC362_portable_v2.zip` supplies prepared data and the tested frozen scientific
snapshot. Use the clone's Python environment to execute the package's runner.
The Gaia data and required external selection code are not in the public clone.

## 1. Check prerequisites

Install Git and 64-bit Python 3.11, then open a new PowerShell window:

```powershell
git --version
py -3.11 --version
```

## 2. Clone the project

```powershell
New-Item -ItemType Directory -Force C:\GCResearch | Out-Null
Set-Location C:\GCResearch
git clone https://github.com/yohan4567/Tidal-Tail-Analysis.git
Set-Location Tidal-Tail-Analysis
```

If this checkout already exists, enter its directory and run `git pull --ff-only`
before starting a new analysis. Retain your working files if Git reports local
changes; do not discard them simply to update.

## 3. Transfer and extract the prepared NGC 362 inputs

Copy `NGC362_portable_v2.zip` from the original PC to this PC using USB, network
transfer or a file-transfer service. The example assumes it is in Downloads:

```powershell
New-Item -ItemType Directory -Force outputs\portable | Out-Null
Expand-Archive -LiteralPath "$HOME\Downloads\NGC362_portable_v2.zip" -DestinationPath outputs\portable
```

Use the version-2 archive. Replace the ZIP path if necessary and use a fresh
extraction destination. Keep the analysis folder in place after starting.

The folder layout should be:

```text
C:\GCResearch\Tidal-Tail-Analysis\
  requirements.txt
  portable_cluster.py
  outputs\portable\NGC362_portable_v2\
    portable_cluster.py
    portable_manifest.json
    requirements.txt
    code\
    inputs\
    external_study\
```

The commands below select the runner inside `outputs\portable\NGC362_portable_v2`.
It locates its manifest, calibration and frozen code beside itself.

## 4. Create the Python environment

Remain in `C:\GCResearch\Tidal-Tail-Analysis`:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r outputs\portable\NGC362_portable_v2\requirements.txt
if ($LASTEXITCODE -ne 0) { throw 'Dependency installation failed.' }
.\.venv\Scripts\python.exe -m pip check
if ($LASTEXITCODE -ne 0) { throw 'Resolve dependency errors before continuing.' }
```

This uses the package's pinned principal dependency versions. Internet is needed
for installation. The prepared input package supplies Gaia, PARSEC and the
completeness map, so those data do not need to be downloaded again.
PowerShell virtual-environment activation is unnecessary.

## 5. Verify the transferred inputs

```powershell
.\.venv\Scripts\python.exe outputs\portable\NGC362_portable_v2\portable_cluster.py verify
if ($LASTEXITCODE -ne 0) { throw 'Input verification failed.' }
```

For this export, the success message is:

```text
Verified 197 files and 15 population nodes for NGC362
```

The package was tested after relocation, including a 16-tracer, 2-Myr integration
and its observation response. These are small software checks; a full grid and
scientific rate validation remain separate.

## 6. Start the analysis

```powershell
.\.venv\Scripts\python.exe -u outputs\portable\NGC362_portable_v2\portable_cluster.py run --out outputs\NGC362_other_pc_v1 --workers 4
```

The results directory must be new on first use. The runner creates it. Keep
PowerShell open and prevent sleep during computation. Four workers is a starting
setting; choose eight or twelve only when the destination PC has enough CPU and
memory headroom. Worker count is not an operating-system RAM limit.

The run uses rates 1, 3, ... 31 solar masses/Myr; lookback times 100, 200, ... 2000
Myr; and three dynamics seeds. It integrates 48 maximum-duration models with
16,000 tracers each. Three escape-MF offsets are processed separately, each with
five constrained population nodes. Compatible dynamics are reused across offsets.

## 7. Check progress

In a second PowerShell window:

```powershell
Set-Location C:\GCResearch\Tidal-Tail-Analysis
.\.venv\Scripts\python.exe outputs\portable\NGC362_portable_v2\portable_cluster.py status --out outputs\NGC362_other_pc_v1
Get-Content outputs\NGC362_other_pc_v1\run_status.json
```

Rerun `status` for an updated snapshot. It reports completed dynamics, full
observation templates and posterior availability. Per-model integration logs
are under `outputs/NGC362_other_pc_v1/offset_0/grid/rate_*/simulation.log`.
Counts are not an estimate of the percentage of time remaining.

## 8. Resume after interruption

After stopping the old process, return to the same checkout and rerun:

```powershell
.\.venv\Scripts\python.exe -u outputs\portable\NGC362_portable_v2\portable_cluster.py run --out outputs\NGC362_other_pc_v1 --workers 4
```

Compatible checkpoints preserve completed work. Shutdown stops execution and
requires a manual restart. Use one writer per result directory. If a partial
observation response fails, inspect its error before retrying; the program
preserves incomplete files rather than silently deleting them.

## 9. Read results

```text
outputs/NGC362_other_pc_v1/
  run_status.json
  offset_0/posterior/posterior_summary.json
  offset_0/posterior/joint_posterior.png
  offset_minus0p5/posterior/
  offset_minus1/posterior/
```

The three offsets are separate escape-population assumptions. The outputs are
conditional inference products; convergence, recovery/coverage and physical-model
checks are required before claiming a validated absolute mass-loss rate.

Updating the Git checkout leaves the package's frozen snapshot unchanged. A new
scientific implementation requires a separately prepared and validated export.
