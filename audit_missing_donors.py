"""Reproduce a failed real response and retain unsupported donor coordinates."""
import argparse,json,os,sys
from pathlib import Path


def run(job,simulation,out):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    j=json.loads(Path(job).read_text(encoding='utf-8-sig'))
    os.environ['GC_TIDAL_STUDY']=j['study'];sys.path.insert(0,j['code'])
    import numpy as np
    import pandas as pd
    from empirical_population import prepare
    from response_model import ErrorLibrary
    from response_group import run_group
    prior=prepare(j['cluster'],j['isochrone'],j['reference'],j['html'],out/'population',offset=0.,order=5)
    original=ErrorLibrary.donor_support
    def diagnostic(self,g,c,x,y):
        features=np.column_stack([g,c,x,y])/self.scale
        distances,indices=self.tree.query(features,k=1)
        bad=distances>self.radius
        if bad.any():
            frame=pd.DataFrame(dict(g=g[bad],g_rp=c[bad],x=x[bad],y=y[bad],nearest_distance=distances[bad]))
            nearest=self.d.iloc[indices[bad]].reset_index(drop=True)
            for col in ('source_id','ra','dec','phot_g_mean_mag','g_rp'):
                frame['nearest_'+col]=nearest[col].to_numpy()
            frame.to_csv(out/'unsupported_pairs.csv',index=False)
            (out/'diagnostic.json').write_text(json.dumps(dict(stage='unsupported_donors_reproduced',
                pairs_in_batch=len(g),unsupported=len(frame),simulation=str(Path(simulation).resolve()),
                min_nearest_distance=float(distances[bad].min()),max_nearest_distance=float(distances[bad].max()),
                feature_scales=self.scale.tolist(),radius=self.radius,physical_rate_measurement_validated=False),indent=2))
        return original(self,g,c,x,y)
    ErrorLibrary.donor_support=diagnostic
    try:
        run_group(Path(j['selection']),Path(simulation),[Path(prior['response_configs'][2])],[out/'response'],study=Path(j['study']))
    except ValueError:
        if (out/'diagnostic.json').exists():return
        raise
    (out/'diagnostic.json').write_text(json.dumps(dict(stage='response_passed_no_missing_donors',physical_rate_measurement_validated=False),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('job','simulation','out'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    run(a.job,a.simulation,a.out)
