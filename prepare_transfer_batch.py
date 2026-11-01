"""Prepare an explicit, versioned cluster transfer without changing live runs."""
import argparse,json,os,shutil,subprocess,sys,time
from pathlib import Path
from datetime import datetime,timezone
from provenance import sha256_file
from selection_bridge import DEFAULT_STUDY

HERE=Path(__file__).resolve().parent
TEN=['NGC6121','NGC6397','NGC6535','NGC6266','NGC7089','NGC6139','NGC6366','NGC6304','NGC6355','NGC6541']
NOTES={
    'NGC6266':'Single-isochrone conditional baseline; known helium populations are not represented. https://doi.org/10.1093/mnras/stu2198',
    'NGC7089':'Single-isochrone conditional baseline; minority metal-rich population and other population structure require sensitivity analysis. https://doi.org/10.1093/mnras/stv2802',
    'NGC6355':'Literature alpha enhancement is recorded but not included in the scaled-solar baseline; differential extinction and CMD agreement require review.',
}

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(d,indent=2),encoding='utf-8');tmp.replace(p)

def copy_field(cid,fields,sources):
    for folder in sources:
        raw=folder/(cid+'_field.csv.gz');meta=folder/(cid+'_field.json')
        if not raw.exists() or not meta.exists():continue
        for suffix in ('_field.csv.gz','_field.json','_missing_tiles.json'):
            p=folder/(cid+suffix);target=fields/p.name
            if p.exists() and not target.exists():shutil.copy2(p,target)
        return True
    return False

def run(root,clusters,previous=None):
    root=Path(root).resolve();study=root/'study';fields=study/'gaia_fields'
    fields.mkdir(parents=True,exist_ok=True);(study/'catalogs').mkdir(exist_ok=True)
    for name in ('step2_query_gaia_fields.py','step3_tail_search.py','skymath.py','catalogs/master_clusters.csv'):
        target=study/name
        if not target.exists():shutil.copy2(DEFAULT_STUDY/name,target)
        if sha256_file(target)!=sha256_file(DEFAULT_STUDY/name):raise ValueError('External source changed: '+name)
    plan=dict(clusters=clusters,analysis_kind='model_comparison',previous_transfer=str(previous) if previous else None,
        excluded_from_previous=['ESO452-SC11','NGC1851','NGC6341','NGC7099'],
        source_sha256=sha256_file(Path(__file__)),physical_rate_measurement_validated=False)
    saved=root/'transfer_plan.json'
    if saved.exists() and read(saved)!=plan:raise ValueError('Changed transfer plan; choose a fresh versioned root')
    dump(saved,plan)
    os.environ['GC_TIDAL_STUDY']=str(study);os.environ['TIDAL_EXTRA_FIELDS']=str(fields)
    states=read(root/'status.json').get('clusters',{}) if (root/'status.json').exists() else {}
    def update(cid,stage,**extra):
        states.setdefault(cid,dict(population_note=NOTES.get(cid),mass_loss_ready=False)).update(stage=stage,
            updated_utc=datetime.now(timezone.utc).isoformat(),**extra)
        dump(root/'status.json',dict(clusters=states,scope='Observed-data preparation; model preparation is tracked separately.'))
        print(cid,stage,extra.get('error',''),flush=True)
    def execute(script,args,log):
        with log.open('a',encoding='utf-8') as f:
            p=subprocess.run([sys.executable,'-u',str(script),*map(str,args)],cwd=HERE,stdout=f,stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
        if p.returncode:raise RuntimeError('\n'.join(log.read_text(encoding='utf-8',errors='replace').splitlines()[-8:]))
    sources=[fields,DEFAULT_STUDY/'gaia_fields']
    if previous:sources.insert(1,Path(previous)/'study/gaia_fields')
    def observed(cid):
        if states.get(cid,{}).get('observation_ready'):return True
        from run_deadline_queue import memory_snapshot
        while True:
            m=memory_snapshot()
            if not m or m['available_gib']>=4:break
            update(cid,'waiting_for_RAM',available_gib=m['available_gib']);time.sleep(30)
        if not copy_field(cid,fields,sources):return False
        dest=root/'observed'/cid;dest.parent.mkdir(exist_ok=True)
        try:
            update(cid,'preparing_observed_selection')
            if not (dest/'REPORT.md').exists():execute(HERE/'run_pipeline.py',['--cluster',cid,'--legacy','--study',study,'--out',dest],root/(cid+'_selection.log'))
            cfg=dest/'adopted_dynamics.json'
            if not cfg.exists():execute(HERE/'make_config.py',[cid,'--catalog',study/'catalogs/master_clusters.csv','--out',cfg],root/(cid+'_config.log'))
            s=read(dest/'selection/selection.json')
            update(cid,'observations_ready_calibration_pending',observation_ready=True,selected_fit_stars=s['n_fit_stars'],
                field_sha256=sha256_file(fields/(cid+'_field.csv.gz')),selection=str(dest/'selection'))
            return True
        except Exception as exc:update(cid,'observation_preparation_failed',error=str(exc),observation_ready=False);return False
    for cid in clusters:update(cid,states.get(cid,{}).get('stage','queued'))
    for cid in clusters:observed(cid)
    # Let an already-running producer finish its own tile downloads and merging.
    # No second writer is started against its files; completed fields are copied.
    if previous:
        finished=Path(previous)/'observed_transfer_manifest.json'
        while not finished.exists():
            for cid in clusters:
                if not states.get(cid,{}).get('observation_ready'):
                    if not observed(cid) and states[cid]['stage']!='observation_preparation_failed':update(cid,'waiting_for_previous_field_downloads')
            time.sleep(30)
        for cid in clusters:observed(cid)
    missing=[c for c in clusters if not (fields/(c+'_field.csv.gz')).exists()]
    if missing:
        for cid in missing:update(cid,'downloading_real_gaia_field')
        try:execute(study/'step2_query_gaia_fields.py',['--all','--release=dr3','--level=6','--tiles=2',*missing],root/'gaia_download.log')
        except Exception as exc:
            for cid in missing:update(cid,'field_download_failed',error=str(exc))
        for cid in missing:
            if not observed(cid):update(cid,'field_missing_or_observation_failed',observation_ready=False)
    dump(root/'observed_transfer_manifest.json',dict(plan=plan,clusters=states,
        finished_utc=datetime.now(timezone.utc).isoformat(),
        files_sha256={str(p.relative_to(root)).replace('\\','/'):sha256_file(p) for folder in (study,root/'observed') for p in folder.rglob('*')
            if p.is_file() and '__pycache__' not in p.parts and not any(part.endswith('.parts') for part in p.parts) and p.suffix not in ('.pyc','.log')}))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--clusters',nargs='+',default=TEN);p.add_argument('--previous-transfer',type=Path);a=p.parse_args()
    if not a.clusters or len(set(a.clusters))!=len(a.clusters):p.error('A nonempty unique cluster list is required')
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    run(a.root,a.clusters,a.previous_transfer.resolve() if a.previous_transfer else None)
