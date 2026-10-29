"""Import one cluster from the project's existing catalogue without changing it."""
import argparse
import csv
import json
from pathlib import Path
from stream_model import HERE, validate
from provenance import sha256_file
from selection_bridge import DEFAULT_STUDY


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cluster')
    p.add_argument('--catalog',type=Path,default=DEFAULT_STUDY/'catalogs/master_clusters.csv')
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    with a.catalog.open(encoding='utf-8-sig',newline='') as f:
        rows=[r for r in csv.DictReader(f) if r['id']==a.cluster]
    if len(rows)!=1: raise ValueError('Expected exactly one matching cluster ID')
    c=json.loads((HERE/'configs/Pal5.json').read_text())
    mapping=dict(ra_deg='ra',dec_deg='dec',distance_kpc='R_sun',pmra_masyr='pmra',
                 pmdec_masyr='pmdec',radial_velocity_kms='rv',mass_msun='mass',half_mass_radius_pc='r_hm_pc')
    c.update({k:float(rows[0][v]) for k,v in mapping.items()})
    c.update(name=a.cluster,source=f'{a.catalog.resolve()}, id={a.cluster}; catalogue snapshot, not fitted stream orbit')
    c['catalogue_provenance']=dict(path=str(a.catalog.resolve()),sha256=sha256_file(a.catalog),columns=mapping,
        original_sources=dict(phase_space='https://people.smp.uq.edu.au/HolgerBaumgardt/globular/orbits_table.txt',
            mass_and_radius='https://people.smp.uq.edu.au/HolgerBaumgardt/globular/combined_table.txt'),
        caveat='Project catalogue snapshot; these parameters are adopted, not independently remeasured or jointly fitted.')
    validate(c)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('x',encoding='utf-8') as f: json.dump(c,f,indent=2)
    print(a.out.resolve())


if __name__=='__main__': main()
