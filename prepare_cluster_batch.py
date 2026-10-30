"""Prepare real-input selections and a provenance audit for the requested sample.

Does not borrow Pal 5 populations or mass conversions for other clusters.
"""
import argparse,json,subprocess,sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
from selection_bridge import DEFAULT_STUDY
from provenance import sha256_file

CLUSTERS=[('M30','NGC7099'),('NGC 288','NGC288'),('Pal 5','Pal5'),('M92','NGC6341'),
          ('NGC 2298','NGC2298'),('M5','NGC5904'),('NGC 5466','NGC5466'),('Omega Cen','NGC5139'),('NGC 3201','NGC3201')]
HERE=Path(__file__).resolve().parent


def prepare(out,workers=2):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    catalog=DEFAULT_STUDY/'catalogs/master_clusters.csv'
    master=pd.read_csv(catalog).set_index('id')
    columns=['ra','dec','R_sun','pmra','pmdec','rv','mass','r_hm_pc','feh','ebv']
    snapshot=master.loc[[i for _,i in CLUSTERS],columns].reset_index()
    snapshot.to_csv(out/'adopted_parameters.csv',index=False)
    def one(item):
        display,cid=item;dest=out/cid;dest.mkdir(exist_ok=True)
        config=dest/'adopted_dynamics.json'
        if not config.exists():
            subprocess.run([sys.executable,'make_config.py',cid,'--catalog',str(catalog),'--out',str(config)],cwd=HERE,check=True,capture_output=True,text=True)
        data=json.loads(config.read_text(encoding='utf-8'))
        source=DEFAULT_STUDY/'gaia_fields'/f'{cid}_field.csv.gz'
        meta=json.loads(source.with_name(f'{cid}_field.json').read_text())
        entry=dict(cluster=cid,display_name=display,source_field=str(source),source_sha256=sha256_file(source),
            field_radius_deg=meta['radius_deg'],field_version='original_v1_explicit',
            dynamics_config=str(config),parameter_catalog_sha256=sha256_file(catalog),
            missing=['Cluster-specific population age/abundance and CMD validation',
                     'Cluster-specific unfiltered Gaia quality/error parent calibration',
                     'Escape-population versus dynamical mass conversion and mass-function assumptions',
                     'Numerical convergence and injection recovery'],
            status='preparing_real_observations',physical_rate_measurement_validated=False)
        try:
            if not (dest/'selection/selection.json').exists():
                with (dest/'selection.log').open('w',encoding='utf-8') as log:
                    result=subprocess.run([sys.executable,'selection_bridge.py','--study',str(DEFAULT_STUDY),
                        '--cluster',cid,'--legacy','--out',str(dest/'selection')],cwd=HERE,stdout=log,stderr=subprocess.STDOUT)
                if result.returncode:raise RuntimeError('Selection failed; see selection.log')
            sel=json.loads((dest/'selection/selection.json').read_text(encoding='utf-8'))
            entry.update(status='observations_ready_calibration_pending',selection=str(dest/'selection'),
                         selection_sha256=sha256_file(dest/'selection/selection.json'),cutflow=sel.get('observation_cutflow'),
                         n_fit_stars=sel['n_fit_stars'],n_fit_bins=sel['n_fit_bins'],cmd_enabled=sel['cmd_enabled'])
            if cid=='Pal5':
                bundle=HERE/'outputs/Pal5_science_calibration/bundle'
                specs=[bundle/f'tilt_{p}_zp_minus0.017.json' for p in ['minus1','plus0','plus1']]
                if all(p.exists() for p in specs):
                    entry.update(status='conditional_grid_inputs_available',response_configs=[str(p) for p in specs],
                        coupled_config=str(HERE/'outputs/Pal5_joint_validation_v1/inputs/dynamics.json'),
                        age_upper_bound_myr=11000.,missing=['Independent convergence and full-grid recovery; population and mass conversion remain conditional'])
            if cid=='NGC5139':entry['additional_limitations']=['Multiple stellar populations require a justified mixture rather than a single Pal 5 isochrone.']
            if cid=='NGC3201':entry['additional_limitations']=['Crowded low-latitude field and differential extinction require explicit validation.']
        except Exception as exc:entry.update(status='failed',error=str(exc))
        (dest/'input_audit.json').write_text(json.dumps(entry,indent=2),encoding='utf-8')
        print(cid,entry['status'],flush=True);return entry
    with ThreadPoolExecutor(max_workers=workers) as pool:entries=list(pool.map(one,CLUSTERS))
    report=dict(schema_version=1,entries=entries,source_sha256=sha256_file(__file__),
        statement='Real Gaia observations and adopted catalogue parameters; no completed formal physical-rate measurements.')
    (out/'batch_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--workers',type=int,default=2)
    a=p.parse_args();prepare(a.out,a.workers)
