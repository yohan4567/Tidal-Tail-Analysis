"""Resume cluster-specific calibration and run remaining conditional survey grids.

No Pal 5 population/error transfer; literature baselines are explicit assumptions.
One cluster's numerical grid runs at a time, with up to two input preparations.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import json
from pathlib import Path
import subprocess
import sys
import threading
import pandas as pd
from provenance import sha256_file
from run_cluster_survey import SAMPLE
from selection_bridge import DEFAULT_STUDY

HERE=Path(__file__).resolve().parent
PRIORITY=['NGC288','NGC6341','NGC7099','NGC5904','NGC2298','NGC5466']


def checked_isochrone(path,age,mh):
    meta=path.with_suffix('.json')
    if not path.exists() or not meta.exists():return False
    m=json.loads(meta.read_text(encoding='utf-8'))
    if m.get('sha256')!=sha256_file(path) or m.get('age_years')!=age or m.get('metallicity_mh')!=mh:
        raise ValueError('Existing isochrone checksum/age/metallicity mismatch')
    return True


def execute(script,args,log):
    log=Path(log);log.parent.mkdir(parents=True,exist_ok=True)
    with log.open('a',encoding='utf-8') as f:
        f.write('\n'+datetime.now(timezone.utc).isoformat()+' '+script+'\n');f.flush()
        result=subprocess.run([sys.executable,'-u',str(HERE/script),*map(str,args)],cwd=HERE,stdout=f,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW if sys.platform=='win32' else 0)
    if result.returncode:
        raise RuntimeError(script+': '+'\n'.join(log.read_text(encoding='utf-8',errors='replace').splitlines()[-6:]))


def run(root,adoptions,workers=10,response_workers=3,prepare_workers=2,clusters=None):
    root=Path(root).resolve();adoptions=Path(adoptions).resolve()
    rules=json.loads(adoptions.read_text(encoding='utf-8'))['clusters']
    catalog=DEFAULT_STUDY/'catalogs/master_clusters.csv';master=pd.read_csv(catalog).set_index('id')
    order=[cid for _,cid in SAMPLE if cid!='Pal5' and (not clusters or cid in clusters)]
    order.sort(key=lambda cid:PRIORITY.index(cid) if cid in PRIORITY else len(PRIORITY))
    lock=threading.RLock();states={};started=datetime.now(timezone.utc).isoformat()
    def update(cid,stage,**extra):
        with lock:
            old=states.setdefault(cid,dict(cluster=cid,physical_rate_measurement_validated=False))
            old.update(stage=stage,updated_utc=datetime.now(timezone.utc).isoformat(),**extra)
            dest=root/cid;dest.mkdir(parents=True,exist_ok=True)
            (dest/'remaining_job.json').write_text(json.dumps(old,indent=2),encoding='utf-8')
            data=dict(started_utc=started,updated_utc=old['updated_utc'],workers=workers,response_workers=response_workers,
                adoption_sha256=sha256_file(adoptions),entries=[states[c] for c in order if c in states])
            temp=root/'remaining_survey_status.tmp';temp.write_text(json.dumps(data,indent=2),encoding='utf-8');temp.replace(root/'remaining_survey_status.json')
            print(cid,stage,flush=True)
    for cid in order:update(cid,'queued')
    def prepare(cid):
        dest=root/cid;selection=dest/'selection';r=rules.get(cid)
        try:
            if not (selection/'selection.json').exists():
                update(cid,'blocked_missing_observed_selection',reason='Real Gaia field/selection not yet available');return None
            if not r or r.get('hold_reason'):
                update(cid,'blocked_population_review',reason=(r or {}).get('hold_reason','No verified literature population age in adopted input table'));return None
            report=json.loads((selection/'selection.json').read_text(encoding='utf-8'))
            params=master.loc[cid];age=r['age_gyr']*1e9;mh=float(params.feh)
            inp=dest/'calibration_inputs';inp.mkdir(exist_ok=True)
            iso=inp/'isochrone.dat'
            if not checked_isochrone(iso,age,mh):
                update(cid,'fetching_official_isochrone',age_gyr=r['age_gyr'],population_source=r['age_source'])
                source=f"{cid}; age: {r['age_source']}; {r['abundance_policy']}"
                execute('fetch_isochrone.py',['--age-years',age,'--mh',mh,'--source',source,'--out',iso],dest/'population_download.log')
            old=HERE/'outputs/nine_cluster_real_inputs_v1'/cid/'parent_random10.csv.gz'
            parent=old if old.exists() and old.with_suffix('.json').exists() else inp/'parent_random10.csv.gz'
            if not parent.exists() or not parent.with_suffix('.json').exists():
                update(cid,'fetching_gaia_calibration_parent')
                for endpoint in ('esa','aip'):
                    try:
                        execute('fetch_parent_tiles.py',['--ra',report['center_ra_deg'],'--dec',report['center_dec_deg'],
                            '--radius',min(10.,report['field_radius_deg']+.05),'--out',parent,'--endpoint',endpoint],dest/'parent_download.log');break
                    except RuntimeError:
                        if endpoint=='aip':raise
            bundle=dest/'science_calibration/bundle'
            if not (bundle/'calibration_manifest.json').exists():
                update(cid,'building_cluster_calibration')
                execute('prepare_science_calibration.py',['--isochrone',iso,'--parent',parent,'--selection',selection,
                    '--ra',report['center_ra_deg'],'--dec',report['center_dec_deg'],'--radius',report['field_radius_deg']+.05,
                    '--ebv',float(params.ebv),'--out',bundle],dest/'calibration.log')
            manifest=json.loads((bundle/'calibration_manifest.json').read_text())
            if manifest['population_source']['sha256']!=sha256_file(iso) or manifest['parent_source']['original_parent_sha256']!=sha256_file(parent):
                raise ValueError('Calibration bundle input checksum mismatch')
            configs=[bundle/f'tilt_{p}_zp_minus0.017.json' for p in ('minus1','plus0','plus1')]
            if not all(p.exists() for p in configs):raise ValueError('Incomplete population sensitivity bundle')
            dynamics=json.loads((dest/'adopted_dynamics.json').read_text())
            dynamics['mass_evolution']=dict(mode='linear_escape',dynamical_to_population_mass=1.,population_mass_basis='isochrone_current_stellar_mass',
                source_description=f'{cid} conditional baseline: assumed f=1 living-star/dynamical escape conversion, fixed Plummer scale; no wind/remnant evolution. Not independently validated.')
            coupled=dest/'coupled_dynamics.json'
            if coupled.exists() and json.loads(coupled.read_text())!=dynamics:raise ValueError('Frozen coupled dynamics changed')
            coupled.write_text(json.dumps(dynamics,indent=2),encoding='utf-8')
            update(cid,'ready_conditional_grid',calibration_bundle=str(bundle),population_assumptions=r,grid=str(dest/'production_grid_shared_v1'))
            return cid,coupled,configs,r['age_gyr']*1000
        except Exception as exc:
            update(cid,'preparation_failed',error=str(exc));return None
    with ThreadPoolExecutor(max_workers=prepare_workers) as pool:
        futures=[pool.submit(prepare,cid) for cid in order]
        for future in as_completed(futures):
            ready=future.result()
            if ready is None:continue
            cid,coupled,specs,age=ready;dest=root/cid;grid=dest/'production_grid_shared_v1'
            try:
                update(cid,'running_conditional_grid')
                execute('production_study.py',['--config',coupled,'--selection',dest/'selection','--response-configs',*specs,
                    '--cluster-age-myr',age,'--workers',workers,'--response-workers',response_workers,
                    '--integrator','cartesian_batch','--batch-size',128,'--out',grid,'--execute'],dest/'production_shared.log')
                assessment=json.loads((grid/'production_assessment.json').read_text())
                update(cid,'completed_conditional_grid_validation_pending',assessment=assessment)
            except Exception as exc:update(cid,'grid_failed',error=str(exc))
    print('Remaining survey queue finished. Read remaining_survey_status.json for completed/blocked/failed stages.',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=HERE/'outputs/cluster_survey_29_v1')
    p.add_argument('--adoptions',type=Path,default=HERE/'configs/survey_populations_v1.json');p.add_argument('--workers',type=int,default=10)
    p.add_argument('--response-workers',type=int,default=3);p.add_argument('--prepare-workers',type=int,default=2);p.add_argument('--clusters',nargs='+')
    a=p.parse_args()
    if min(a.workers,a.response_workers,a.prepare_workers)<1:p.error('Worker counts must be positive')
    run(a.root,a.adoptions,a.workers,a.response_workers,a.prepare_workers,a.clusters)
