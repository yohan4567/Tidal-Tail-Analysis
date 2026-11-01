"""Exercise relocated runtime paths with an existing small real-input simulation."""
import argparse,json,os
from pathlib import Path
import portable_cluster as portable

def run(bundle,simulation,out):
    bundle,simulation,out=map(lambda p:Path(p).resolve(),(bundle,simulation,out))
    info=portable.verify(bundle)
    portable.prepare(bundle,out/'bound_run',info)
    # Scientific imports occur only after verify selects the package's snapshot.
    from response_group import run_group
    from validate_science import check_response
    prior=portable.read(out/'bound_run/offset_0/population_prior.json')
    response=out/'relocated_response'
    run_group(bundle/'inputs/selection',simulation,[prior['response_configs'][2]],[response],study=bundle/'external_study')
    meta=check_response(response)
    portable.dump(out/'validation.json',dict(cluster=info['cluster'],stage='relocated_small_response_passed',
        bundle_manifest_sha256=portable.digest(bundle/'portable_manifest.json'),template_sha256=meta['template_sha256'],
        validator_sha256=portable.digest(Path(__file__)),physical_rate_measurement_validated=False))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('bundle','simulation','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
    run(a.bundle,a.simulation,a.out)
