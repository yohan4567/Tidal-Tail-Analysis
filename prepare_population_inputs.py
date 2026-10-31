"""Fetch missing cluster-specific isochrones; preserve existing/frozen inputs."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd
from fetch_isochrone import fetch
from population_catalog import save,HERE


def run(root):
    root=Path(root).resolve();catalog=save(root)
    master=pd.read_csv(root/'adopted_parameters.csv').set_index('id');entries=[]
    def update(entry):
        entry['updated_utc']=datetime.now(timezone.utc).isoformat()
        (root/'population_preparation_status.json').write_text(json.dumps(dict(entries=entries),indent=2))
        print(entry['cluster'],entry['status'],entry.get('reason',''),flush=True)
    for record in catalog['clusters']:
        cid=record['cluster'];age=record['age'];iso=(HERE/record['isochrone']['path']).resolve()
        entry=dict(cluster=cid,status='checking');entries.append(entry)
        try:
            if age.get('hold_reason') or age.get('age_gyr') is None:
                entry.update(status='blocked_multi_population',reason=age.get('hold_reason','No single age'));update(entry);continue
            if iso.exists() or iso.with_suffix('.json').exists():
                entry.update(status='preserved_existing_input',reason=record['isochrone']['status']);update(entry);continue
            mh=age.get('measured_mh',float(master.loc[cid,'feh']))
            if not np.isfinite(mh):raise ValueError('No explicit adopted metallicity')
            entry.update(status='fetching_official_isochrone',age_gyr=age['age_gyr'],requested_mh=mh);update(entry)
            fetch(age['age_gyr']*1e9,mh,iso,source=f"{cid}; {age['age_source']}; age status: {age.get('age_status','existing literature adoption')}; scaled-solar baseline, no inferred escape MF")
            entry['status']='isochrone_ready';update(entry)
        except Exception as exc:
            entry.update(status='blocked_or_failed',reason=str(exc));update(entry)
    save(root)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=HERE/'outputs/cluster_survey_29_v1');a=p.parse_args();run(a.root)
