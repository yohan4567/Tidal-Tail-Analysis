"""Versioned NGC6366 calibration and M80 donor-support repairs with real tests."""
import argparse,json,os,shutil,subprocess,sys,time,zipfile
from pathlib import Path
from datetime import datetime,timezone
from provenance import sha256_file
from repair_seven_clusters import worker_lock
from transfer_status_io import write_json

HERE=Path(__file__).resolve().parent

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))

def dump(p,data):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);write_json(p,data,status_file=p.name=='status.json')

def run(root,run_after=None):
    root=Path(root).resolve();states={}
    def status(cid,stage,**extra):
        d=dict(stage=stage,updated_utc=datetime.now(timezone.utc).isoformat(),physical_rate_measurement_validated=False)
        d.update(extra);states.setdefault(cid,{}).update(d);dump(root/'status.json',dict(clusters=states))
        print(cid,stage,extra.get('error',''),flush=True)
    def memory(cid):
        from run_deadline_queue import memory_snapshot
        while (m:=memory_snapshot()) and m['available_gib']<8:
            status(cid,'waiting_for_RAM',available_gib=m['available_gib']);time.sleep(30)
    def command(script,args,cid):
        memory(cid)
        log=root/(cid+'.log')
        with log.open('a',encoding='utf-8') as stream:
            r=subprocess.run([sys.executable,'-u',str(script),*map(str,args)],cwd=HERE,stdout=stream,stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
        if r.returncode:raise RuntimeError('\n'.join(log.read_text(encoding='utf-8',errors='replace').splitlines()[-10:]))
    from recover_compact_repairs import stage
    from queue_repaired_seven import make_job
    ready={}
    for cid in ('NGC6366','NGC6093'):
        try:
            dest=root/cid
            if cid=='NGC6366':
                status(cid,'rebuilding_training_partition_and_exact_discrete_screen')
                stage(HERE/'outputs/netr2/NGC6366/prepared',dest,HERE/'outputs/rc3',pool=True,pool_screening='exact_beta_binomial_v2')
            else:
                original=HERE/'outputs/s7r2/NGC6093';settings=read(original/'bright_repair.json')
                source=Path(settings['source'])/cid;selection=read(source/'selection/selection.json')
                os.environ['GC_TIDAL_STUDY']=settings['study']
                extra=root/'inputs/M80_G12_18_census.csv.gz'
                status(cid,'fetching_complete_G12_18_error_donor_census')
                from repair_parent_support import fetch_parent
                fetch_parent(extra,selection,magnitude_range=(12.,18.),sampling_modulus=1,prefer_tiles=True,
                    tile_level=7,tiles_per_request=1,before_merge=lambda:memory(cid))
                memory(cid);stage(original,dest,HERE/'outputs/rc3',pool=False)
                manifest=dest/'donor_extension.json'
                if not manifest.exists():
                    import pandas as pd
                    from selection_bridge import load_study
                    from prepare_science_calibration import spatial_parent
                    newsettings=read(dest/'bright_repair.json');bundle=Path(newsettings['reference']).parent
                    oldspec=read(newsettings['reference']);oldfile=Path(oldspec['errors']['path'])
                    a=pd.read_csv(oldfile,low_memory=False);b=pd.read_csv(extra,low_memory=False)
                    T=load_study(Path(settings['study']),selection['legacy']);b=spatial_parent(b,T,selection)
                    if not ((b.phot_g_mean_mag>=12)&(b.phot_g_mean_mag<18)).all():raise ValueError('Supplement outside declared G range')
                    combined=pd.concat([a,b],ignore_index=True).drop_duplicates('source_id',keep='first')
                    path=bundle/'parent_errors_extended.csv.gz';combined.to_csv(path,index=False,compression='gzip')
                    meta=read(oldfile.with_suffix('.json'));meta.update(n_sources=len(combined),sha256=sha256_file(path),
                        source_description=meta['source_description']+'; error-library-only full Gaia DR3 G=12..18 census supplement, with identical sky/core/neighbour exclusions; quality map unchanged',
                        sampling_strata=[dict(g_range=[0,18],fraction=1.),dict(g_range=[18,21],fraction=.1)],
                        error_support_extension=dict(original_parent_sha256=sha256_file(oldfile),supplement_metadata=read(extra.with_suffix('.json')),
                            new_unique_sources=len(combined)-len(a),quality_counts_unchanged=True,donor_metric_and_radius_unchanged=True,
                            caveat='G-conditioned finite-neighbour empirical donor distribution; the richer sampling changes that distribution and still requires transfer/systematic validation.'))
                    dump(path.with_suffix('.json'),meta)
                    for cfg in bundle.glob('tilt_*.json'):
                        s=read(cfg);s['errors'].update(path=str(path),metadata=str(path.with_suffix('.json')));dump(cfg,s)
                    dump(manifest,meta['error_support_extension'])
                    del a,b,combined
            if not (dest/'repair_manifest.json').exists():
                status(cid,'testing_real_response');command(HERE/'recover_compact_repairs.py',['--finish',dest],cid)
            jobpath=root/cid/'local/job.json';make_job(dest,cid,jobpath)
            if cid=='NGC6093':
                from prepare_exact_local_retries import physical_config
                job=read(jobpath);original_cfg=HERE/'outputs/cluster_survey_29_v1/NGC6093/coupled_dynamics.json'
                if physical_config(read(job['config']))!=physical_config(read(original_cfg)):raise ValueError('M80 dynamics no longer physically identical')
                job['config']=str(original_cfg);job['frozen_sha256'][str(original_cfg)]=sha256_file(original_cfg)
                grid=HERE/'outputs/s7r2/NGC6093/local_v3/results/offset_0/grid';job['reuse_grid']=str(grid)
                retryjob=jobpath.with_name('reuse_job.json')
                if retryjob.exists() and read(retryjob)!=job:raise ValueError('M80 reuse job changed')
                dump(retryjob,job);jobpath=retryjob
                status(cid,'auditing_all_48_existing_simulations')
                command(HERE/'audit_grid_support.py',['--job',jobpath,'--grid',grid,'--out',root/cid/'grid_support'],cid)
                audit=read(root/cid/'grid_support/audit.json')
                if not audit['complete'] or audit['checked']!=48 or any(not r['passed'] for r in audit['results']):raise ValueError('M80 exhaustive support check incomplete')
                # Also rerun the exact previously failing case with all native
                # observation/calibration draws, not just the support-only audit.
                status(cid,'testing_previously_failed_full_duration_response')
                cfg=read(root/cid/'grid_support/population/population_prior.json')['response_configs'][2]
                command(Path(job['code'])/'response_group.py',['--selection',job['selection'],'--run',grid/'rate_27_seed_3037/simulation',
                    '--configs',cfg,'--outputs',root/cid/'failed_case_response','--study',job['study']],cid)
            package=root/'pkgs'/(cid+'_v3');validation=root/'v'/cid
            if not package.with_suffix('.zip').exists():
                status(cid,'exporting_package');command(HERE/'portable_cluster.py',['export','--repair',dest,'--cluster',cid,
                    '--out',package,'--analysis-kind','model_comparison'],cid)
            if not (validation/'validation.json').exists():
                status(cid,'checking_relocated_response');command(HERE/'validate_transfer_package.py',
                    ['--bundle',package,'--simulation',dest/'preflight/simulation','--out',validation],cid)
            if read(validation/'validation.json')['stage']!='relocated_small_response_passed':raise ValueError('Relocation validation missing')
            archive=package.with_suffix('.zip')
            with zipfile.ZipFile(archive) as z:
                if z.testzip() is not None:raise ValueError('ZIP CRC verification failed')
            entry=dict(archive=str(archive),archive_sha256=sha256_file(archive),bytes=archive.stat().st_size,zip_crc_verified=True,
                       validation=str(validation/'validation.json'),job=str(jobpath),analysis_kind='model_comparison',physical_rate_measurement_validated=False)
            ready[cid]=entry;dump(root/'READY_PACKAGES.json',dict(ready_packages=ready));status(cid,'conditional_model_package_ready',**entry)
            if cid=='NGC6366':
                p=HERE/'outputs/other_pc_all_remaining_20261101_v1/READY_PACKAGES.json';index=read(p)
                index['ready_packages'][cid]=entry;index['updated_utc']=datetime.now(timezone.utc).isoformat();dump(p,index)
        except Exception as exc:status(cid,'repair_failed',error=str(exc))
    dump(root/'preparation_finished.json',dict(clusters=states))
    if run_after is not None and 'NGC6093' in ready:
        cid='NGC6093';entry=ready[cid]
        def cancelled():
            if not (root/cid/'local_retry_cancelled.json').exists():return False
            status(cid,'local_rerun_cancelled_by_user',**entry)
            return True
        if cancelled():return
        while read(Path(run_after)/'status.json').get('stage')!='queue_finished_check_individual_statuses':
            if cancelled():return
            status(cid,'ready_waiting_for_existing_queue',**entry);time.sleep(30)
        if cancelled():return
        try:
            status(cid,'running_repaired_model_comparisons',**entry)
            command(HERE/'run_priority_model.py',['--job',entry['job'],'--out',root/cid/'results','--workers',12],cid)
            if read(root/cid/'results/run_status.json')['stage']!='completed_model_comparisons_validation_pending':raise ValueError('Missing completed marker')
            status(cid,'completed_model_comparisons_validation_pending',**entry)
        except Exception as exc:status(cid,'analysis_failed',error=str(exc))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--run-after',type=Path);a=p.parse_args()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    with worker_lock(a.root):run(a.root,a.run_after)
