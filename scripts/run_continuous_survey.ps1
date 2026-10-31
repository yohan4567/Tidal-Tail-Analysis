param(
    [int]$Workers = 12,
    [double]$MemoryBudgetGiB = 24,
    [switch]$RetryFailed
)
$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $repo
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$queueArgs = @('-u', 'run_deadline_queue.py', '--workers', $Workers, '--memory-budget-gib', $MemoryBudgetGiB)
if ($RetryFailed) { $queueArgs += '--retry-failed' }
& (Join-Path $repo '.venv/Scripts/python.exe') @queueArgs
exit $LASTEXITCODE
