param(
    [Parameter(Mandatory=$true)][int]$WaitPid,
    [Parameter(Mandatory=$true)][string]$OutputRoot,
    [int]$Workers=12
)
$ErrorActionPreference='Stop'
$repo=Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $repo
New-Item -ItemType Directory -Force $OutputRoot | Out-Null
$root=(Resolve-Path -LiteralPath $OutputRoot).Path
function Save-State([string]$Stage) {
    [ordered]@{stage=$Stage;pid=$PID;waiting_for_pid=$WaitPid;workers=$Workers;
        order=@('NGC6144','NGC5466','NGC5904','NGC6093','NGC2298','NGC3201');
        updated_utc=[DateTime]::UtcNow.ToString('o')} | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $root 'scheduled_run.json') -Encoding UTF8
}
Save-State 'waiting_for_current_M107'
$existing=Get-Process -Id $WaitPid -ErrorAction SilentlyContinue
if ($existing) {
    $started=$existing.StartTime
    while ($true) {
        $current=Get-Process -Id $WaitPid -ErrorAction SilentlyContinue
        if (-not $current -or $current.StartTime -ne $started) { break }
        Start-Sleep -Seconds 30
    }
}
$env:OPENBLAS_NUM_THREADS='1';$env:OMP_NUM_THREADS='1';$env:MKL_NUM_THREADS='1'
Save-State 'running_priority_six'
$ErrorActionPreference='Continue'
& (Join-Path $repo '.venv/Scripts/python.exe') -u (Join-Path $repo 'priority_queue.py') --root $root --workers $Workers *> (Join-Path $root 'queue.log')
$code=$LASTEXITCODE
$ErrorActionPreference='Stop'
if ($code -ne 0) { Save-State 'controller_failed';exit $code }
Save-State 'priority_pass_finished_check_cluster_statuses'
$repair=Join-Path $repo 'outputs/input_repairs_20261101_v1'
# Preserve other work after the six requested priorities. The preparation process
# may still be downloading inputs; wait without launching a second numerical grid.
Save-State 'waiting_for_remaining_input_preparation'
$preparation=Get-Process -Id 17752 -ErrorAction SilentlyContinue
if ($preparation) {
    $preparationStart=$preparation.StartTime
    while ($true) {
        $current=Get-Process -Id 17752 -ErrorAction SilentlyContinue
        if (-not $current -or $current.StartTime -ne $preparationStart) { break }
        Start-Sleep -Seconds 30
    }
}
Save-State 'continuing_other_prepared_repairs'
$ErrorActionPreference='Continue'
& (Join-Path $repo '.venv/Scripts/python.exe') -u (Join-Path $repo 'dispatch_repair_queue.py') --repair $repair --workers $Workers *> (Join-Path $root 'remaining_repairs.log')
$ErrorActionPreference='Stop'
Save-State 'continuing_remaining_original_targets'
$ErrorActionPreference='Continue'
& (Join-Path $repo '.venv/Scripts/python.exe') -u (Join-Path $repo 'run_deadline_queue.py') --workers $Workers --memory-budget-gib 24 *> (Join-Path $root 'remaining_targets.log')
$ErrorActionPreference='Stop'
Save-State 'all_passes_finished_check_cluster_statuses'
