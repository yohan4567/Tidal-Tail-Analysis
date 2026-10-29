"""Compare supplied, predeclared stream templates on the identical observed count grid."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from fit_rate import analyse, analyse_ensemble,require_free_amplitude
from provenance import sha256_file


def compare(paths,out):
    if len(paths)<2:
        raise ValueError('At least two predeclared templates required')
    reference=None; reference_description=None
    results=[]
    for path in map(Path,paths):
        with np.load(path,allow_pickle=False) as d:
            n,s,b=[d[k] for k in ('counts','signal','background')]
            ref=float(d['reference_mdot']); desc=str(d['selection_description'].item())
            draws=d['signal_draws'] if 'signal_draws' in d else None
            synthetic=bool(d['synthetic_calibration'].item()) if 'synthetic_calibration' in d else False
        try:
            description=json.loads(desc)
        except json.JSONDecodeError:
            raise ValueError('Template comparison requires structured selection provenance')
        if not isinstance(description,dict) or 'selection_sha256' not in description:
            raise ValueError('Selection identity required')
        require_free_amplitude(description)
        identity={k:description.get(k) for k in ('selection_sha256','observation_identity','mass_basis')}
        if not identity.get('observation_identity'):
            raise ValueError('Regenerate response templates with observation identity before comparing models')
        if reference is None:
            reference=(n.copy(),b.copy()); reference_description=identity
        elif not np.array_equal(n,reference[0]) or not np.array_equal(b,reference[1]) or identity!=reference_description:
            raise ValueError('Models must use identical observed bins, background and selection/mass definitions')
        result=analyse(n,s,b,ref) if draws is None else analyse_ensemble(n,draws,b,ref)
        result.update(template=str(path.resolve()),template_sha256=sha256_file(path),synthetic_calibration=synthetic)
        results.append(result)
    minimum=min(r['negative_log_likelihood'] for r in results)
    for result in results:
        result['delta_twice_nll']=2*(result['negative_log_likelihood']-minimum)
    results.sort(key=lambda r:r['negative_log_likelihood'])
    out=Path(out); out.mkdir(parents=True,exist_ok=False)
    summary=[{k:r[k] for k in ('template','amplitude','background_scale','mdot_msun_myr','negative_log_likelihood','delta_twice_nll','synthetic_calibration')} for r in results]
    pd.DataFrame(summary).to_csv(out/'comparison.csv',index=False)
    report=dict(models=results,observation_identity=reference_description,
                interpretation='Relative conditional profile scores on one common dataset; not Bayes factors or automatic detection significance.',
                limitations=['Fixed candidate models and common background shape; searched model/setting trials require separate calibration.',
                             'Dynamic parameters and unrepresented calibration systematics are not integrated automatically.'])
    (out/'comparison.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--templates',type=Path,nargs='+',required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    r=compare(a.templates,a.out)
    print(json.dumps([dict(template=x['template'],delta_twice_nll=x['delta_twice_nll']) for x in r['models']],indent=2))
