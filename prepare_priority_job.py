"""Prepare one isolated repair without changing live survey scientific modules."""
import argparse,json,os,shutil,subprocess,sys
from pathlib import Path
import pandas as pd
from provenance import sha256_file
from selection_bridge import DEFAULT_STUDY
from selection_model import QualityMap,parent_query,download_parent
from quality_repair import validate_pooling
from validate_science import absolute_spec

HERE=Path(__file__).resolve().parent


def dump(path,data):Path(path).write_text(json.dumps(data,indent=2),encoding='utf-8')
def prepare(cluster,out):
    out=Path(out).resolve();job=out/'job.json'
    if job.exists():return job
    if out.exists():raise ValueError('Incomplete preparation retained; inspect it and choose a fresh version')
    out.mkdir(parents=True);source=HERE/'outputs/cluster_survey_29_v1'/cluster
    parent_repair=HERE/'outputs/input_repairs_20261101_v1';base=json.loads((parent_repair/'repair_manifest.json').read_text())
    for file,sha in base['frozen_sha256'].items():
        if sha256_file(file)!=sha:raise ValueError('Base repair source changed: '+file)
    reference=source/'science_calibration/bundle/tilt_plus0_zp_minus0.017.json'
    if cluster in ('NGC5904','NGC3201'):
        selection=json.loads((source/'selection/selection.json').read_text());bright=out/'bright_parent.csv.gz'
        q=parent_query(selection['center_ra_deg'],selection['center_dec_deg'],selection['field_radius_deg']+.05,0.,12.,1)
        download_parent(q,bright,100000)
        from bright_calibration import prepare as extend
        extend(parent_repair,cluster,bright,out/'bright_extension')
        info=json.loads((out/'bright_extension/bright_repair.json').read_text());reference=Path(info['reference']);base_code=Path(info['snapshot'])
    else:base_code=Path(base['snapshot'])
    code=out/'code';shutil.copytree(base_code,code,ignore=shutil.ignore_patterns('__pycache__','outputs'))
    # Explicit loader and hashes for conditional spatial pooling in this snapshot.
    p=code/'selection_model.py';text=p.read_text(encoding='utf-8');needle='        return cls(arrays,meta,draw_seed)'
    assert text.count(needle)==1
    text=text.replace(needle,"        if meta.get('quality_pooling'):\n            from quality_repair import SpatialPooledQualityMap\n            return SpatialPooledQualityMap(arrays,meta,draw_seed)\n"+needle)
    p.write_text(text,encoding='utf-8')
    for filename in ('response_model.py','run_joint_study.py'):
        p=code/filename;text=p.read_text(encoding='utf-8')
        text=text.replace("'selection_model.py',","'selection_model.py','quality_repair.py',",1);p.write_text(text,encoding='utf-8')
    for filename in ('likelihood_comparison.py','run_priority_model.py'):shutil.copy2(HERE/filename,code/filename)
    spec=absolute_spec(reference);quality=QualityMap.load(spec['quality']);validation=None
    if (quality.arrays['n_available']<quality.minimum).any():
        held=Path(reference).parent/'quality_holdout_cells.csv';candidate,validation=validate_pooling(quality,pd.read_csv(held))
        if not validation['screening_passed']:raise ValueError('Pooled quality failed held-out screen; additional calibration data required')
        candidate.save(out/'quality');spec.update(quality=str(out/'quality'),quality_repair=candidate.metadata['quality_pooling'])
        dump(out/'pooling_validation.json',validation)
    # Use the official full-coverage coarse map for the known multi-map holes.
    if cluster=='NGC5466':
        map_path=HERE/'outputs/cache/gaia/allsky_M10_hpx7.hdf5'
        if not map_path.exists():raise ValueError('Official hpx7 completeness map required')
        spec['completeness'].update(path=str(map_path),mode='hpx7',
            source_description='Official GaiaUnlimited TCG hpx7 all-sky map; ICRS query convention matched to upstream implementation; separate coarse-resolution repair of multi-map gaps.')
    spec['cache_dir']=str(out/'response_cache');reference=out/'reference_response.json';dump(reference,spec)
    grid=source/'empirical_direct_v1/offset_0/grid';reuse=grid if (grid/'plan.json').exists() else None
    rules=json.loads((HERE/'configs/survey_populations_v1.json').read_text())['clusters']
    for cid,extra in json.loads((HERE/'configs/population_supplements_v1.json').read_text())['clusters'].items():rules[cid]=dict(rules.get(cid,{}),**extra)
    info=dict(cluster=cluster,code=str(code),study=str(DEFAULT_STUDY.resolve()),config=str(source/'coupled_dynamics.json'),
        selection=str(source/'selection'),isochrone=str(source/'calibration_inputs/isochrone.dat'),
        reference=str(reference),html=str(HERE/'outputs/empirical_mf_source/parameter.html'),reuse_grid=str(reuse) if reuse else None,
        cluster_age_myr=rules[cluster]['age_gyr']*1000,validation=validation,
        target_parameter_posterior_computed=False,physical_rate_measurement_validated=False)
    files=[*code.glob('*.py'),*list((code/'configs').rglob('*.json')),reference,Path(info['config']),Path(info['isochrone']),
        Path(info['isochrone']).with_suffix('.json'),Path(info['html']),Path(info['selection'])/'selection.json',Path(info['selection'])/'maps.npz',
        Path(spec['quality'])/'quality.json',Path(spec['quality'])/'quality.npz',Path(spec['errors']['path']),Path(spec['errors']['metadata']),Path(spec['completeness']['path'])]
    info['frozen_sha256']={str(p.resolve()):sha256_file(p) for p in files};dump(job,info);return job


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cluster',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    print(prepare(a.cluster,a.out))
