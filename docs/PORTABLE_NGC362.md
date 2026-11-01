# Run the prepared NGC 362 analysis on another Windows PC

Use the locally supplied `NGC362_portable_v2.zip`. This archive contains the
prepared observations, actual Gaia calibration inputs, official completeness
map, corrected weighted populations, required external selection code and a
frozen analysis-code snapshot. It contains no virtual environment, old numerical
grid or mass-loss result. The raw data and external code in this private transfer
package are not uploaded to the public GitHub repository.

The package is about 93 MiB. Transfer it to the destination computer by USB,
local network or your preferred file-transfer service. Do not copy a virtual
environment between computers.

## 1. Install Python and extract

Install 64-bit Python 3.11. In PowerShell, confirm:

```powershell
py -3.11 --version
```

Assuming the ZIP is in Downloads:

```powershell
New-Item -ItemType Directory -Force C:\GCResearch | Out-Null
Expand-Archive -LiteralPath "$HOME\Downloads\NGC362_portable_v2.zip" -DestinationPath C:\GCResearch
Set-Location C:\GCResearch\NGC362_portable_v2
```

If the archive is elsewhere, replace its path. Use a new extraction directory;
keep this run outside a cloud-synchronised folder. No Git clone is required for
this prepared package.

## 2. Create the environment

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
if ($LASTEXITCODE -ne 0) { throw 'Dependency installation failed.' }
.\.venv\Scripts\python.exe -m pip check
if ($LASTEXITCODE -ne 0) { throw 'Resolve dependency errors before continuing.' }
```

The package pins the principal dependencies to the source environment's
versions. Internet access is needed to install them. The prepared analysis does
not need to redownload Gaia, PARSEC or the completeness map.

## 3. Verify the transferred inputs

```powershell
.\.venv\Scripts\python.exe portable_cluster.py verify
if ($LASTEXITCODE -ne 0) { throw 'Input verification failed; do not start the grid.' }
```

Expected message: `Verified ... files and 15 population nodes for NGC362`.
This checks file hashes, common photometric support, positive weights and unit
represented stellar mass. It does not establish physical validity of the rate.

The runner binds current file locations at runtime. Historical source paths in
unchanged provenance records may still name the original computer; do not edit
them or the immutable input tables to replace user names.

## 4. Start the analysis

```powershell
.\.venv\Scripts\python.exe -u portable_cluster.py run --workers 4
```

Four workers is a conservative starting setting, not a measured RAM guarantee.
On a computer with sufficient cores and 32 GiB or more RAM, you can choose eight
workers; use twelve only with sufficient CPU and memory headroom. Change the
number in the command, not the scientific grid or frozen input files.

The command verifies inputs, creates `results/portable_binding.json`, prepares
local response paths and starts the conditional grid. It uses:

- Escape-rate assumptions: 1, 3, ... 31 solar masses per Myr (16 values).
- Lookback times: 100, 200, ... 2000 Myr (20 values).
- Three dynamics seeds: 48 integrations to the maximum lookback time, with
  16,000 tracers each. Smaller lookback times use release-time subsets, with an
  automatic independent-window numerical check.
- Three separate escape-MF offsets: 0, -0.5, -1, with five constrained global-MF
  nodes per offset. Later offsets reuse verified dynamics, not new integrations.
- 32 observation draws and eight shared calibration states.

The nominal rate/time grid has 320 locations per population/seed scenario.
These are model hypotheses, not 320 independent observational measurements.
No full-run completion time is promised by the small portability test.

## 5. Read progress

Open another PowerShell window:

```powershell
Set-Location C:\GCResearch\NGC362_portable_v2
.\.venv\Scripts\python.exe portable_cluster.py status
Get-Content results\run_status.json
```

Run `status` whenever you want a fresh snapshot. It reports completed dynamics,
full observation templates and posterior availability for each offset. Individual
integration logs are in `results/offset_0/grid/rate_*/simulation.log`.
These counts are not a percentage of elapsed time or scientific validation.

## 6. Resume and locate results

After stopping an interrupted run, rerun the same command from the same package
location with the same output path. Never run two writers in the same results
directory. A shutdown stops the program; it does not automatically restart.
Compatible checkpoints preserve completed work. A persistent failed partial
response must be inspected before resuming; the runner does not silently delete
it. Keep the package in place after binding; choose a fresh `--out` if moving it.

```powershell
.\.venv\Scripts\python.exe -u portable_cluster.py run --workers 4
```

Successful conditional posterior products are under:

```text
results/
  offset_0/posterior/
  offset_minus0p5/posterior/
  offset_minus1/posterior/
```

Look for `posterior_summary.json` and `joint_posterior.png`. The offsets remain
separate assumptions about the escaping population. Results still require
response-support, convergence, recovery/coverage and physical-applicability
validation before being presented as an absolute physical mass-loss rate.

## Population assumptions preserved in this transfer

The weighted populations use the supplied PARSEC isochrone (adopted age 10.75
Gyr and scaled-solar [M/H] = -1.26), the recorded external global-MF constraint
and the three escape-MF offsets. The abundance mapping is a baseline
approximation. Normalization is to represented living-star current mass;
unrepresented remnants, winds and masses outside the supplied isochrone are not
automatically recovered by running the grid. The dynamical conversion remains
the documented assumption f = 1.

The export was checked after extraction to a different directory containing a
space. Verification, path binding and a 16-tracer, 2-Myr integration completed.
This small test checks portability; it does not replace the full experiment.

The first export's observation test found model giants at G about 11.8, outside
the original G=12..21 calibration. Version 2 extends calibration with a complete
Gaia bright-star census below G=12, keeping the independent 10% faint sample and
placing the sampling-design boundary at a quality-bin edge. The original
selection and photometric nodes are retained. Use the version-2 archive; the
version-1 export did not pass the observation-response preflight.
