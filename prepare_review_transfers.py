"""Export observations and explicit input-review tasks for the nine held targets."""
import argparse,json,os,shutil,subprocess,sys,time,zipfile
from pathlib import Path
from datetime import datetime,timezone
from importlib.metadata import version
from provenance import sha256_file

HERE=Path(__file__).resolve().parent
REVIEWS={
    'ESO452-SC11':'Global MF is indirect, without a measured mass interval. Reconcile missing catalogue abundance/reddening with the literature; do not label a model-based MF as measured.',
    'NGC1851':'Review a multi-population composition and weight model; a single isochrone has not been validated for absolute-rate interpretation.',
    'NGC6341':'Requested MH=-2.31, returned MH=-2.19174. Adopt a supported, composition-consistent isochrone and preserve its actual returned parameters.',
    'NGC7099':'Requested MH=-2.27, returned MH=-2.19174. Adopt a supported, composition-consistent isochrone and preserve its actual returned parameters.',
    'NGC5139':'Requires a multi-age, multi-composition population with appropriate weights; the mean age is not a replacement for that distribution.',
    'NGC5286':'Review the iron/s-process populations and their photometric weights; single-population conversion is not validated.',
    'NGC6273':'Review the metallicity distribution and population weights; single-population conversion is not validated.',
    'NGC6656':'Review the complex populations and abundance-model dependence before choosing isochrones and mixture weights.',
    'NGC7078':'Existing isochrone table has MH=-2.19174 and lacks its download metadata. Recover a verified, composition-consistent model and provenance; do not create replacement metadata for this unverified table.',
}

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,d):Path(p).write_text(json.dumps(d,indent=2),encoding='utf-8')

def run(root,wait_for,previous):
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True)
    state=root/'review_transfer_status.json'
    def update(stage,**extra):
        dump(state,dict(stage=stage,updated_utc=datetime.now(timezone.utc).isoformat(),clusters=list(REVIEWS),**extra))
        print(stage,flush=True)
    update('waiting_for_ten_cluster_observation_preparation')
    while not (Path(wait_for)/'observed_transfer_manifest.json').exists():time.sleep(30)
    update('preparing_remaining_observations')
    with (root/'observations.log').open('a',encoding='utf-8') as log:
        result=subprocess.run([sys.executable,'-u',str(HERE/'prepare_transfer_batch.py'),'--root',str(root),
            '--previous-transfer',str(Path(previous).resolve()),'--clusters',*REVIEWS],cwd=HERE,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
    if result.returncode:raise RuntimeError('Observation producer failed; inspect observations.log')
    update('packaging_review_inputs')
    package=root/'review_work_package';package.mkdir(exist_ok=False)
    code=package/'code';code.mkdir()
    for p in HERE.glob('*.py'):shutil.copy2(p,code/p.name)
    shutil.copytree(HERE/'configs',code/'configs');shutil.copy2(HERE/'LICENSE',package/'LICENSE')
    study=package/'study';(study/'catalogs').mkdir(parents=True);(study/'gaia_fields').mkdir()
    for name in ('step2_query_gaia_fields.py','step3_tail_search.py','skymath.py','catalogs/master_clusters.csv'):
        shutil.copy2(root/'study'/name,study/name)
    for cid in REVIEWS:
        for suffix in ('_field.csv.gz','_field.json','_missing_tiles.json'):
            p=root/'study/gaia_fields'/(cid+suffix)
            if p.exists():shutil.copy2(p,study/'gaia_fields'/p.name)
    observations=read(root/'status.json')['clusters'];status={}
    for cid,reason in REVIEWS.items():
        dest=package/'inputs'/cid;dest.mkdir(parents=True)
        source=HERE/'outputs/cluster_survey_29_v1'/cid/'calibration_inputs'
        for suffix in ('dat','json'):
            p=source/('isochrone.'+suffix)
            if p.exists():shutil.copy2(p,dest/p.name)
        obs=root/'observed'/cid
        if obs.exists():shutil.copytree(obs,dest/'observed',ignore=shutil.ignore_patterns('__pycache__','*.log'))
        status[cid]=dict(observation_ready=observations.get(cid,{}).get('observation_ready',False),
            conditional_mass_loss_run_ready=False,physical_rate_measurement_validated=False,required_review=reason)
    shutil.copy2(HERE/'outputs/input_repairs_20261101_v1/mf_source.html',package/'mf_source.html')
    dump(package/'READINESS.json',status)
    (package/'requirements.txt').write_text('\n'.join(f'{p}=={version(p)}' for p in
        ('galpy','astropy','numpy','scipy','matplotlib','pandas','astropy-healpix','h5py','requests'))+'\n',encoding='utf-8')
    (package/'RUN_GUIDE.md').write_text('# Observation and input-review work package\n\n'
        'These nine targets are assigned to the other PC. Read READINESS.json before interpreting results.\n'
        'This package supports observation/tail analysis and review of the supplied inputs; it is not a ready absolute-rate grid.\n'
        'Known mismatched or unverifiable isochrones are retained as evidence and must not be used as verified models.\n\n'
        'Install requirements.txt in a new Python 3.11 environment. To regenerate one observation analysis from the package root:\n\n'
        '```powershell\n.\\.venv\\Scripts\\python.exe code/run_pipeline.py --cluster NGC1851 --legacy --study study --out results/NGC1851\n```\n\n'
        'Choose a fresh output directory; inspect failed targets separately. Model corrections require scientific review, not removal of support or provenance checks.\n',encoding='utf-8')
    files={str(p.relative_to(package)).replace('\\','/'):sha256_file(p) for p in package.rglob('*') if p.is_file()}
    dump(package/'files_manifest.json',dict(files_sha256=files,source_sha256=sha256_file(Path(__file__))))
    archive=package.with_suffix('.zip')
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=1) as z:
        for p in package.rglob('*'):
            if p.is_file():z.write(p,str(Path(package.name)/p.relative_to(package)))
    with zipfile.ZipFile(archive) as z:
        bad=z.testzip()
        if bad:raise ValueError('Archive failed CRC verification: '+bad)
    update('observation_and_review_package_ready',archive=str(archive),archive_sha256=sha256_file(archive),
        readiness=status,physical_rate_measurement_validated=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','wait-for','previous-transfer'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
    run(a.root,a.wait_for,a.previous_transfer)
