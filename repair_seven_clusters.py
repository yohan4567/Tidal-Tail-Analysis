"""Repair seven requested GC inputs in isolated versions, then exercise real responses."""
import argparse,json,os,shutil,subprocess,sys,time,zipfile
from pathlib import Path
from datetime import datetime,timezone
from provenance import sha256_file
from transfer_status_io import write_json
from repair_parent_support import fetch_parent
from repair_four_populations import stage_code,read,cmd_review
from fetch_isochrone import validate_returned_model
from population_model import read_isochrone

HERE=Path(__file__).resolve().parent
CLUSTERS=('NGC5904','NGC6093','NGC6341','NGC7099','NGC7078','NGC6266','NGC6139')

def dump(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    write_json(p,d,status_file=p.name=='status.json')

def run(root):
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True)
    composition=root/'composition';cal=root/'calibration'
    states=read(root/'status.json').get('clusters',{}) if (root/'status.json').exists() else {}
    def update(cid,stage,**extra):
        if stage!='repair_failed':states.setdefault(cid,{}).pop('error',None)
        states.setdefault(cid,{}).update(stage=stage,updated_utc=datetime.now(timezone.utc).isoformat(),
            physical_rate_measurement_validated=False,**extra)
        dump(root/'status.json',dict(clusters=states));print(cid,stage,extra.get('error',''),flush=True)
    def memory(cid,minimum=8):
        from run_deadline_queue import memory_snapshot
        while True:
            m=memory_snapshot()
            if not m or m['available_gib']>=minimum:return
            update(cid,'waiting_for_RAM',available_gib=m['available_gib'],required_gib=minimum)
            time.sleep(30)
    def command(script,args,cid):
        with (root/(cid+'_model.log')).open('a',encoding='utf-8') as f:
            r=subprocess.run([sys.executable,'-u',str(script),*map(str,args)],cwd=HERE,stdout=f,stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
        if r.returncode:raise RuntimeError('\n'.join((root/(cid+'_model.log')).read_text(encoding='utf-8',errors='replace').splitlines()[-12:]))
    for cid in CLUSTERS:
        if states.get(cid,{}).get('stage')!='conditional_model_package_ready':update(cid,'queued')
    from selection_bridge import DEFAULT_STUDY
    study=cal/'study';(study/'catalogs').mkdir(parents=True,exist_ok=True)
    for name in ('step3_tail_search.py','skymath.py','catalogs/master_clusters.csv'):
        dest=study/name
        if not dest.exists():shutil.copy2(DEFAULT_STUDY/name,dest)
    os.environ['GC_TIDAL_STUDY']=str(study)
    base,code=stage_code(cal,study)
    source_roots={}
    for cid in CLUSTERS:
        try:
            if cid in ('NGC6341','NGC7099','NGC7078'):
                if not (composition/'cases'/cid/'review.json').exists():raise ValueError('Reviewed composition inputs are not yet complete')
                source_root=composition
            else:
                source_root=cal
                source=HERE/'outputs/cluster_survey_29_v1'/cid
                observed=source if cid in ('NGC5904','NGC6093') else HERE/'outputs/other_pc_ten_20261101_v1/observed'/cid
                dest=cal/'observed'/cid;dest.mkdir(parents=True,exist_ok=True)
                if not (dest/'selection').exists():shutil.copytree(observed/'selection',dest/'selection')
                if not (dest/'adopted_dynamics.json').exists():shutil.copy2(observed/'adopted_dynamics.json',dest/'adopted_dynamics.json')
                inp=dest/'calibration_inputs';inp.mkdir(exist_ok=True)
                for suffix in ('dat','json'):
                    p=inp/('isochrone.'+suffix)
                    if not p.exists():shutil.copy2(source/'calibration_inputs'/p.name,p)
                cfg=read(dest/'adopted_dynamics.json')
                cfg['mass_evolution']=dict(mode='linear_escape',dynamical_to_population_mass=1.,
                    population_mass_basis='isochrone_current_stellar_mass',
                    source_description='Conditional represented living-star baseline f=1; no remnant/wind correction; fixed Plummer scale.')
                dump(dest/'coupled_dynamics.json',cfg)
                s=read(dest/'selection/selection.json')
                if sha256_file(study/'step3_tail_search.py')!=s['selection_file_sha256']:raise ValueError('Frozen selector differs from staged study')
                if sha256_file(s['source_field'])!=s['source_field_sha256']:raise ValueError('Observed field checksum mismatch')
                meta=read(inp/'isochrone.json')
                validate_returned_model(read_isochrone(inp/'isochrone.dat'),meta['age_years'],meta['metallicity_mh'])
                if not (inp/'population_adoption.json').exists():
                    dump(inp/'population_adoption.json',dict(cluster=cid,source_isochrone=meta,
                        scientific_scope='Existing scaled-solar conditional population retained; this support repair does not establish an escape MF or validate total absolute mass loss.',
                        population_caveat='Known helium populations are not represented by the single-isochrone baseline.' if cid=='NGC6266' else 'Composition, extinction and escape-MF assumptions require sensitivity checks.'))
            source_roots[cid]=source_root
        except Exception as exc:update(cid,'repair_failed',error=str(exc))
    # Downloading precedes memory-intensive calibrations. Complete components and
    # disjoint tile caches survive restart; selected science fields are never parents.
    parents_ready=[]
    for cid in CLUSTERS:
        if cid not in source_roots:continue
        if states.get(cid,{}).get('stage')=='conditional_model_package_ready':continue
        try:
            dest=source_roots[cid]/'observed'/cid;inp=dest/'calibration_inputs'
            s=read(dest/'selection/selection.json')
            candidates=[];reuse=[]
            for old in [HERE/'outputs/nine_cluster_real_inputs_v1'/cid,
                        HERE/'outputs/input_repairs_20261101_v1'/cid/'calibration_inputs',
                        HERE/'outputs/other_pc_ten_20261101_v1/observed'/cid/'calibration_inputs',
                        HERE/'outputs/cluster_survey_29_v1'/cid/'calibration_inputs']:
                candidates.append(old/'parent_random10.csv.gz')
                reuse.extend([old,old/'parent_random10.csv.gz.parts'])
            update(cid,'verifying_or_downloading_faint_parent')
            faint=fetch_parent(inp/'parent_random10.csv.gz',s,candidates=candidates,reuse_dirs=reuse,
                before_merge=lambda:memory(cid))
            update(cid,'verifying_or_downloading_bright_census',faint_parent_stars=faint['n_sources'])
            bright=fetch_parent(inp/'bright_parent.csv.gz',s,bright=True,
                candidates=[p.parent/'bright_parent.csv.gz' for p in candidates])
            update(cid,'uncut_calibration_parents_ready',faint_parent_stars=faint['n_sources'],bright_parent_stars=bright['n_sources'])
            parents_ready.append(cid)
        except Exception as exc:update(cid,'repair_failed',error=str(exc))
    for cid in parents_ready:
        try:
            source_root=source_roots[cid];dest=source_root/'observed'/cid;inp=dest/'calibration_inputs'
            # Every component is isolated; completed source manifests are untouched.
            source_base=source_root/'model_base'
            folder=root/'model_repairs'/cid;folder.mkdir(parents=True,exist_ok=True)
            attempts=sorted(folder.glob('attempt_*'))
            extension=next((p for p in reversed(attempts) if (p/'bright_repair.json').exists()),None)
            if extension is None:
                extension=folder/f'attempt_{len(attempts)+1:03d}'
                memory(cid);update(cid,'building_repaired_calibration')
                command(code/'bright_calibration.py',['--repair',source_base,'--cluster',cid,
                    '--bright',inp/'bright_parent.csv.gz','--out',extension],cid)
            if not (extension/'preflight.json').exists():
                for p in inp.iterdir():
                    if p.name in ('population_adoption.json','cmd_review.json','cmd_review.png') or p.name.startswith('isochrone_alpha'):
                        shutil.copy2(p,extension/cid/'calibration_inputs'/p.name)
                memory(cid);update(cid,'testing_small_real_input_response')
                command(Path(read(extension/'bright_repair.json')['snapshot'])/'bright_calibration.py',
                        ['--finalize','--out',extension],cid)
            state=read(extension/'input_status.json')[cid]
            j=dict(cluster=cid,code=read(extension/'repair_manifest.json')['snapshot'],
                study=str(source_root/'study') if source_root==composition else str(study),
                config=str(dest/'coupled_dynamics.json'),selection=str(dest/'selection'),
                isochrone=str(inp/'isochrone.dat'),reference=state['reference'],html=str(extension/'mf_source.html'),
                reuse_grid=None,cluster_age_myr=read(inp/'isochrone.json')['age_years']/1e6,
                target_parameter_posterior_computed=False,physical_rate_measurement_validated=False,
                frozen_sha256=read(extension/'repair_manifest.json')['frozen_sha256'])
            j['frozen_sha256'].update({str((dest/'selection/maps.npz').resolve()):sha256_file(dest/'selection/maps.npz')})
            if cid in ('NGC5904','NGC6093'):
                old=HERE/'outputs/cluster_survey_29_v1'/cid/'empirical_direct_v1/offset_0/grid'
                if (old/'plan.json').exists():j['reuse_grid']=str(old)
            job=folder/'job.json'
            if not job.exists():dump(job,j)
            elif read(job)!=j:raise ValueError('Existing job identity differs; choose a new repair version')
            # In M80 the original failure occurred with the full-duration tail.
            # A tiny smoke test alone cannot demonstrate that it has been repaired.
            if cid in ('NGC5904','NGC6093'):
                if not j['reuse_grid']:raise ValueError('Existing full-duration dynamics required for support validation')
                memory(cid);update(cid,'testing_existing_full_duration_response')
                command(HERE/'run_priority_model.py',['--job',job,'--out',folder/'support_validation','--workers',1,'--validate-only'],cid)
            package=root/'packages'/(cid+'_portable_v1')
            if not package.with_suffix('.zip').exists():
                memory(cid);update(cid,'exporting_conditional_package')
                command(HERE/'portable_cluster.py',['export','--repair',extension,'--cluster',cid,'--out',package,'--analysis-kind','model_comparison'],cid)
            validation=root/'package_validation'/cid
            if not (validation/'validation.json').exists():
                memory(cid);update(cid,'checking_relocated_package')
                command(HERE/'validate_transfer_package.py',['--bundle',package,'--simulation',extension/'preflight/simulation','--out',validation],cid)
            v=read(validation/'validation.json')
            if v['stage']!='relocated_small_response_passed':raise ValueError('Relocation response has not passed')
            archive=package.with_suffix('.zip')
            with zipfile.ZipFile(archive) as z:
                bad=z.testzip()
                if bad is not None:raise ValueError('ZIP CRC failure: '+bad)
            index=read(root/'READY_PACKAGES.json') if (root/'READY_PACKAGES.json').exists() else dict(ready_packages={})
            entry=dict(archive=str(archive),bytes=archive.stat().st_size,archive_sha256=sha256_file(archive),
                validation=str(validation/'validation.json'),zip_crc_verified=True,analysis_kind='model_comparison',
                physical_rate_measurement_validated=False,job=str(job))
            index['ready_packages'][cid]=entry;dump(root/'READY_PACKAGES.json',index)
            update(cid,'conditional_model_package_ready',**entry)
        except Exception as exc:update(cid,'repair_failed',error=str(exc))
    dump(root/'finished.json',dict(clusters=states,finished_utc=datetime.now(timezone.utc).isoformat()))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    run(a.root)
