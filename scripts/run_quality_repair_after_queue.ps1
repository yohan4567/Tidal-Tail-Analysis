param(
    [Parameter(Mandatory=$true)][string]$RepairPath,
    [Parameter(Mandatory=$true)][int]$QueuePid,
    [int]$Workers=12
)
$ErrorActionPreference='Stop'
if ($Workers -lt 1 -or $Workers -gt 12) { throw 'Use 1–12 workers.' }
$repo=Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $repo
$repair=(Resolve-Path -LiteralPath $RepairPath).Path
$statePath=Join-Path $repair 'scheduled_run.json'
function Write-RepairState([string]$Stage,[int]$Code=0) {
    [ordered]@{stage=$Stage;pid=$PID;waiting_for_pid=$QueuePid;workers=$Workers;
        updated_utc=[DateTime]::UtcNow.ToString('o');exit_code=$Code;
        physical_rate_measurement_validated=$false} | ConvertTo-Json | Set-Content -LiteralPath $statePath -Encoding UTF8
}
Write-RepairState 'waiting_for_existing_queue'
$existing=Get-Process -Id $QueuePid -ErrorAction SilentlyContinue
if ($existing) {
    $started=$existing.StartTime
    while ($true) {
        $current=Get-Process -Id $QueuePid -ErrorAction SilentlyContinue
        if (-not $current -or $current.StartTime -ne $started) { break }
        Start-Sleep -Seconds 30
    }
}
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
Write-RepairState 'running_repaired_grid'
try {
    # Native stderr can contain harmless scientific-library warnings. The
    # Python exit code, rather than a PowerShell NativeCommandError, decides.
    $ErrorActionPreference='Continue'
    & (Join-Path $repo '.venv/Scripts/python.exe') -u (Join-Path $repair 'frozen_code/repair_quality_calibration.py') run --repair $repair --workers $Workers *> (Join-Path $repair 'full_run.log')
    $result=$LASTEXITCODE
    $ErrorActionPreference='Stop'
    if ($result -ne 0) { Write-RepairState 'failed' $result;exit $result }
    Write-RepairState 'completed_conditional_posteriors_validation_pending'
} catch {
    $_ | Out-String | Add-Content -LiteralPath (Join-Path $repair 'scheduler_error.log')
    Write-RepairState 'failed' 1
    throw
}
