param(
    [Parameter(Mandatory=$true)][string]$RepairPath,
    [Parameter(Mandatory=$true)][int[]]$WaitPids,
    [int]$Workers=12
)
$ErrorActionPreference='Stop'
if ($Workers -lt 1 -or $Workers -gt 12) { throw 'Use 1–12 workers.' }
$repo=Split-Path -Parent $PSScriptRoot
$repair=(Resolve-Path -LiteralPath $RepairPath).Path
$dispatcher=Join-Path $repair 'dispatch/dispatch_repair_queue.py'
if (-not (Test-Path -LiteralPath $dispatcher)) {
    New-Item -ItemType Directory -Force (Split-Path -Parent $dispatcher) | Out-Null
    Copy-Item -LiteralPath (Join-Path $repo 'dispatch_repair_queue.py') -Destination $dispatcher
    Copy-Item -LiteralPath (Join-Path $repo 'provenance.py') -Destination (Join-Path (Split-Path -Parent $dispatcher) 'provenance.py')
}
$excluded=@()
$overrides=Join-Path $repair 'dispatch_overrides.json'
if (Test-Path -LiteralPath $overrides) { $excluded=@((Get-Content -LiteralPath $overrides -Raw | ConvertFrom-Json).excluded_clusters) }
$targets=@()
foreach ($targetPid in $WaitPids) {
    $process=Get-Process -Id $targetPid -ErrorAction SilentlyContinue
    if ($process) { $targets+=@{id=$targetPid;started=$process.StartTime} }
}
function Save-State([string]$Stage,[int]$Code=0) {
    [ordered]@{stage=$Stage;pid=$PID;workers=$Workers;waiting_for_pids=@($targets | ForEach-Object {$_.id});excluded_clusters=$excluded;
        updated_utc=[DateTime]::UtcNow.ToString('o');exit_code=$Code;
        physical_rate_measurement_validated=$false} | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $repair 'scheduled_run.json') -Encoding UTF8
}
Save-State 'waiting_for_existing_jobs_and_input_preparation'
while ($targets.Count -gt 0) {
    $remaining=@()
    foreach ($target in $targets) {
        $process=Get-Process -Id $target.id -ErrorAction SilentlyContinue
        if ($process -and $process.StartTime -eq $target.started) { $remaining+=$target }
    }
    $targets=$remaining
    Save-State 'waiting_for_existing_jobs_and_input_preparation'
    if ($targets.Count -gt 0) { Start-Sleep -Seconds 30 }
}
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
Save-State 'running_repaired_clusters'
Set-Location -LiteralPath $repo
$ErrorActionPreference='Continue'
& (Join-Path $repo '.venv/Scripts/python.exe') -u $dispatcher --repair $repair --workers $Workers *> (Join-Path $repair 'inference.log')
$result=$LASTEXITCODE
$ErrorActionPreference='Stop'
if ($result -ne 0) { Save-State 'runner_failed' $result;exit $result }
# The Python runner deliberately continues after individual cluster failures.
# Its per-cluster inference_status.json files determine scientific completion.
Save-State 'finished_check_per_cluster_results'
