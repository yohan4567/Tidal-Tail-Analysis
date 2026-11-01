"""Prepare real observed fields for the requested nine-cluster PC transfer.

Population holds remain explicit; observed-ready never means mass-loss-ready.
"""
import argparse,json,os,shutil,subprocess,sys,zipfile
from pathlib import Path
from datetime import datetime,timezone
from provenance import sha256_file
from selection_bridge import DEFAULT_STUDY

HERE=Path(__file__).resolve().parent
CLUSTERS=['ESO452-SC11','NGC1851','NGC6121','NGC6139','NGC6341','NGC6366','NGC6397','NGC7099','NGC6535']
HOLDS={'ESO452-SC11':'Indirect MF only; missing catalogue composition/reddening must be sourced. No direct-MF absolute-rate run.',
    'NGC1851':'Multiple-population composition/weight model requires review.',
    'NGC6341':'Requested and returned PARSEC metallicity mismatch; replacement model required.',
    'NGC7099':'Requested and returned PARSEC metallicity mismatch; replacement model required.'}


def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,d):Path(p).write_text(json.dumps(d,indent=2),encoding='utf-8')


def run(root):
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True);study=root/'study';study.mkdir(exist_ok=True)
    fields=study/'gaia_fields';fields.mkdir(exist_ok=True);catalog=study/'catalogs';catalog.mkdir(exist_ok=True)
    for name in ('step2_query_gaia_fields.py','step3_tail_search.py','skymath.py'):
        target=study/name
        if not target.exists():shutil.copy2(DEFAULT_STUDY/name,target)
        if sha256_file(target)!=sha256_file(DEFAULT_STUDY/name):raise ValueError('External code changed during transfer preparation')
    if not (catalog/'master_clusters.csv').exists():shutil.copy2(DEFAULT_STUDY/'catalogs/master_clusters.csv',catalog/'master_clusters.csv')
    for cid in CLUSTERS:
        for suffix in ('_field.csv.gz','_field.json','_missing_tiles.json'):
            p=DEFAULT_STUDY/'gaia_fields'/(cid+suffix);target=fields/p.name
            if p.exists() and not target.exists():shutil.copy2(p,target)
    os.environ['GC_TIDAL_STUDY']=str(study);os.environ['TIDAL_EXTRA_FIELDS']=str(fields)
    states=read(root/'status.json').get('clusters',{}) if (root/'status.json').exists() else {}
    def update(cid,stage,**extra):
        states.setdefault(cid,dict(population_hold=HOLDS.get(cid),mass_loss_ready=False)).update(stage=stage,updated_utc=datetime.now(timezone.utc).isoformat(),**extra)
        dump(root/'status.json',dict(clusters=states,scope='Real observations and sources; model/calibration preparation tracked separately.'))
        print(cid,stage,extra.get('error',''),flush=True)
    def execute(script,args,log):
        with log.open('a',encoding='utf-8') as f:
            f.write('\n'+datetime.now(timezone.utc).isoformat()+'\n');f.flush()
            r=subprocess.run([sys.executable,'-u',str(script),*map(str,args)],cwd=HERE,stdout=f,stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
        if r.returncode:raise RuntimeError('\n'.join(log.read_text(encoding='utf-8',errors='replace').splitlines()[-7:]))
    for cid in CLUSTERS:update(cid,states.get(cid,{}).get('stage','queued'))
    def observed(cid):
        dest=root/'observed'/cid;dest.parent.mkdir(exist_ok=True)
        if not (fields/(cid+'_field.csv.gz')).exists() or not (fields/(cid+'_field.json')).exists():return False
        try:
            update(cid,'preparing_observed_selection')
            if not (dest/'REPORT.md').exists():
                execute(HERE/'run_pipeline.py',['--cluster',cid,'--legacy','--study',study,'--out',dest],root/(cid+'_selection.log'))
            config=dest/'adopted_dynamics.json'
            if not config.exists():execute(HERE/'make_config.py',[cid,'--catalog',catalog/'master_clusters.csv','--out',config],root/(cid+'_config.log'))
            selection=read(dest/'selection/selection.json')
            update(cid,'observations_ready_model_review_pending' if cid in HOLDS else 'observations_ready_calibration_pending',
                observation_ready=True,selected_fit_stars=selection['n_fit_stars'],selection=str(dest/'selection'),
                field_sha256=sha256_file(fields/(cid+'_field.csv.gz')))
            return True
        except Exception as exc:update(cid,'observation_preparation_failed',error=str(exc));return False
    missing=[]
    for cid in CLUSTERS:
        if not (fields/(cid+'_field.csv.gz')).exists():missing.append(cid);update(cid,'waiting_for_real_gaia_field')
        else:observed(cid)
    if missing:
        for cid in missing:update(cid,'downloading_real_gaia_field')
        # The existing downloader queries only sky/G/PM windows; its outputs are
        # observations, never the uncut quality-calibration parent.
        try:execute(study/'step2_query_gaia_fields.py',['--all','--release=dr3','--level=6','--tiles=2',*missing],root/'gaia_download.log')
        except Exception as exc:
            for cid in missing:update(cid,'field_download_failed',error=str(exc))
        for cid in missing:
            if not observed(cid):update(cid,'field_missing_or_observation_failed',observation_ready=False)
    # Package the read-only observed stage separately, preserving each cluster's
    # ID, source field, query metadata, selection and unresolved model holds.
    manifest=dict(created_utc=datetime.now(timezone.utc).isoformat(),clusters=states,
        observed_stage_only=True,mass_loss_ready=False,
        preparation_source_sha256=sha256_file(Path(__file__)),
        files_sha256={str(p.relative_to(root)).replace('\\','/'):sha256_file(p) for folder in (study,root/'observed') if folder.exists() for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ('.log','.pyc')})
    dump(root/'observed_transfer_manifest.json',manifest)
    print('Observed input preparation pass finished. Consult status.json; mass-loss input validation is a separate stage.',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    run(a.root)
