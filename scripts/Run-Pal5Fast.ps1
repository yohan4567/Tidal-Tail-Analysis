param([switch]$PlanOnly)
$projectRoot=Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $projectRoot
$env:MPLCONFIGDIR=Join-Path $projectRoot 'outputs/.mplcache'
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:NUMEXPR_NUM_THREADS='1'
$pal5Arguments=@('-u','production_study.py',
 '--config','outputs/Pal5_joint_validation_v1/inputs/dynamics.json',
 '--selection','outputs/cluster_survey_29_v1/Pal5/selection',
 '--response-configs',
 'outputs/Pal5_science_calibration/bundle/tilt_minus1_zp_minus0.017.json',
 'outputs/Pal5_science_calibration/bundle/tilt_plus0_zp_minus0.017.json',
 'outputs/Pal5_science_calibration/bundle/tilt_plus1_zp_minus0.017.json',
 '--cluster-age-myr','11000','--workers','10','--response-workers','2',
 '--integrator','cartesian_batch','--batch-size','128',
 '--rate-min','1','--rate-max','31','--rate-step','2',
 '--reuse-grid','outputs/cluster_survey_29_v1/Pal5/production_grid_fast_v1',
 '--out','outputs/cluster_survey_29_v1/Pal5/production_grid_odd_v1')
if (-not $PlanOnly) { $pal5Arguments+='--execute' }
& (Join-Path $projectRoot '.venv/Scripts/python.exe') @pal5Arguments
exit $LASTEXITCODE
