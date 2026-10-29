"""Import isochrones and normalize stellar counts to an explicitly defined mass basis."""
import argparse
from dataclasses import dataclass
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd

from provenance import sha256_file

MASS_BASES=('initial_formed_mass','isochrone_current_stellar_mass')
MAG_COLUMNS=('M_G','M_BP','M_RP')


def read_isochrone(path):
    """Read canonical CSV or PARSEC's commented-header whitespace format."""
    path=Path(path)
    text=path.read_text(encoding='utf-8-sig')
    if path.suffix.lower()=='.csv':
        return pd.read_csv(io.StringIO(text))
    header=None
    for line in text.splitlines():
        if line.lstrip().startswith('#') and 'Mini' in line.split() and 'Gmag' in line.split():
            header=line.lstrip()[1:].split()
    if header is None:
        raise ValueError('No PARSEC Mini/Gmag header found. Supply a Gaia-photometry isochrone or canonical CSV.')
    table=pd.read_csv(io.StringIO(text),sep=r'\s+',comment='#',names=header)
    if table.empty:
        raise ValueError('Empty isochrone')
    return table


def build_population(table,*,mass_basis,weight_mode='int_imf',alpha=None,
                     source_description,age_years,metallicity_dex,synthetic=False,
                     escape_mass_tilt=0.):
    if mass_basis not in MASS_BASES or not str(source_description).strip():
        raise ValueError('Explicit supported mass basis and source description required')
    if not np.isfinite([age_years,metallicity_dex]).all() or age_years<=0:
        raise ValueError('Positive age and finite metallicity required')
    # A file containing several ages/compositions cannot silently become one population.
    for col in ('logAge','MH','Zini'):
        if col in table and table[col].nunique(dropna=False)!=1:
            raise ValueError(f'Multiple {col} values; split the isochrone file first')
    if 'logAge' in table and not np.isclose(float(table.logAge.iloc[0]),np.log10(age_years),atol=1e-5,rtol=0):
        raise ValueError('Declared age disagrees with isochrone logAge')
    if 'MH' in table and not np.isclose(float(table.MH.iloc[0]),metallicity_dex,atol=1e-4,rtol=0):
        raise ValueError('Declared metallicity disagrees with isochrone MH')
    aliases={'Mini':'initial_mass_msun','Mass':'current_mass_msun','Gmag':'M_G',
             'G_BPmag':'M_BP','G_RPmag':'M_RP','label':'phase'}
    t=table.rename(columns={k:v for k,v in aliases.items() if v not in table}).copy()
    # CMD may append a dark remnant/termination point, not a complete remnant
    # population. Including its wide initial-mass interval biases normalization.
    excluded=0
    if 'phase' in t:
        keep=pd.to_numeric(t.phase,errors='coerce')!=9
        excluded+=int((~keep).sum()); t=t.loc[keep].copy()
    if 'logL' in t:
        keep=t.logL>-9
        excluded+=int((~keep).sum()); t=t.loc[keep].copy()
    required=['initial_mass_msun','current_mass_msun',*MAG_COLUMNS]
    if any(k not in t for k in required) or len(t)<2:
        raise ValueError('Need at least two isochrone nodes with initial/current masses and G/BP/RP absolute magnitudes')
    if not np.isfinite(t[required].to_numpy(float)).all() or (t.initial_mass_msun<=0).any() or (t.current_mass_msun<=0).any():
        raise ValueError('Nonfinite isochrone or nonpositive masses')
    t=t.sort_values('initial_mass_msun',kind='stable').reset_index(drop=True)
    m=t.initial_mass_msun.to_numpy(float)
    if weight_mode=='int_imf':
        if 'int_IMF' not in t or not np.isfinite(t.int_IMF).all():
            raise ValueError('int_IMF column required; do not assign equal node probabilities')
        cumulative=t.int_IMF.to_numpy(float)
        delta=np.diff(cumulative)
        if (delta<0).any() or delta.sum()<=0:
            raise ValueError('Cumulative IMF must be nondecreasing and nonconstant')
        # Trapezoidal node weights integrate over the supplied initial-mass interval.
        weights=.5*(np.r_[0.,delta]+np.r_[delta,0.])
        initial_normalizer=1.  # PARSEC int_IMF is already per initially formed Msun.
    elif weight_mode=='power_law_initial_coordinate':
        if alpha is None or not np.isfinite(alpha) or (np.diff(m)<=0).any():
            raise ValueError('Finite alpha and strictly increasing initial mass needed for power-law quadrature')
        density=m**(-alpha)
        widths=.5*(np.r_[0.,np.diff(m)]+np.r_[np.diff(m),0.])
        weights=density*widths
        initial_normalizer=float(np.sum(weights*m))
    else:
        raise ValueError('Unsupported weight mode')
    if not np.isfinite(escape_mass_tilt):
        raise ValueError('Finite escape-mass tilt required')
    if escape_mass_tilt!=0:
        if mass_basis!='isochrone_current_stellar_mass':
            raise ValueError('Escape-mass tilt requires current stellar mass normalization')
        # A declared sensitivity family, not a measurement of the escaping mass function.
        # Retain isochrone evolutionary phase weights instead of treating each row equally.
        weights*=np.power(t.current_mass_msun.to_numpy(float)/.5,escape_mass_tilt)
    denom=initial_normalizer if mass_basis=='initial_formed_mass' else float(np.sum(weights*t.current_mass_msun))
    if not np.isfinite(denom) or denom<=0:
        raise ValueError('Invalid mass normalization')
    t['number_per_msun']=weights/denom
    t['bp_rp']=t.M_BP-t.M_RP
    t['g_rp']=t.M_G-t.M_RP
    if 'phase' not in t:
        t['phase']='unspecified'
    extinction_columns=[k for k in ('A_G','A_BP','A_RP') if k in t]
    if extinction_columns and len(extinction_columns)!=3:
        raise ValueError('Supply all three per-node extinctions or none')
    if extinction_columns and (not np.isfinite(t[extinction_columns].to_numpy(float)).all() or (t[extinction_columns]<0).any().any()):
        raise ValueError('Invalid per-node extinction')
    t=t.loc[t.number_per_msun>0,[*required,'phase','bp_rp','g_rp','number_per_msun',*extinction_columns]].reset_index(drop=True)
    meta=dict(schema_version=2,mass_basis=mass_basis,source_description=source_description,
              age_years=float(age_years),metallicity_dex=float(metallicity_dex),weight_mode=weight_mode,alpha=alpha,
              escape_mass_tilt=float(escape_mass_tilt),
              excluded_dark_termination_nodes=excluded,
              escape_mass_function='Baseline isochrone weights times (current_mass / 0.5 Msun)**escape_mass_tilt, renormalized to represented living stellar mass',
              synthetic_calibration=bool(synthetic),n_nodes=len(t),
              represented_initial_mass_range=[float(m.min()),float(m.max())],
              represented_current_mass_per_basis_msun=float(np.sum(t.current_mass_msun*t.number_per_msun)),
              stars_per_basis_msun=float(t.number_per_msun.sum()),
              assumptions=['Distribution is the supplied isochrone plus the declared initial-coordinate weights; not automatically a measured escape mass function.',
                           'No automatic remnant, binary or missing low-mass extension; mass basis covers the declared model only.',
                           'int_IMF quadrature spans only supplied isochrone nodes, with no extrapolation.'])
    return Population(t,meta)


@dataclass
class Population:
    table: pd.DataFrame
    metadata: dict

    def __post_init__(self):
        cols=[*MAG_COLUMNS,'number_per_msun']
        if any(k not in self.table for k in cols) or self.table.empty:
            raise ValueError('Population needs M_G, M_BP, M_RP and number_per_msun')
        if not np.isfinite(self.table[cols].to_numpy(float)).all() or (self.table.number_per_msun<0).any() or self.table.number_per_msun.sum()<=0:
            raise ValueError('Invalid population magnitudes/weights')
        if self.metadata.get('mass_basis') not in MASS_BASES:
            raise ValueError('Population mass_basis missing or unsupported')
        if self.metadata.get('mass_basis')=='isochrone_current_stellar_mass':
            if 'current_mass_msun' not in self.table or not np.isfinite(self.table.current_mass_msun).all() or (self.table.current_mass_msun<=0).any():
                raise ValueError('Current-stellar-mass basis needs current_mass_msun for normalization verification')
            if not np.isclose(np.sum(self.table.current_mass_msun*self.table.number_per_msun),1.,rtol=1e-6,atol=1e-8):
                raise ValueError('Current stellar masses times number_per_msun must sum to one Msun')
        if not self.metadata.get('source_description','').strip():
            raise ValueError('Population source description required')
        extinction=[k for k in ('A_G','A_BP','A_RP') if k in self.table]
        if extinction and (len(extinction)!=3 or not np.isfinite(self.table[extinction].to_numpy(float)).all() or (self.table[extinction]<0).any().any()):
            raise ValueError('Population per-node extinction must have all three finite nonnegative bands')

    def save(self,path):
        path=Path(path)
        path.mkdir(parents=True,exist_ok=False)
        self.table.to_csv(path/'population.csv',index=False)
        meta=dict(self.metadata,table_sha256=sha256_file(path/'population.csv'))
        (path/'population.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')

    @classmethod
    def load(cls,path):
        path=Path(path)
        meta=json.loads((path/'population.json').read_text(encoding='utf-8'))
        if meta.get('table_sha256')!=sha256_file(path/'population.csv'):
            raise ValueError('Population table hash mismatch; rebuild metadata instead of silently editing weights')
        return cls(pd.read_csv(path/'population.csv'),meta)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--isochrone',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--mass-basis',choices=MASS_BASES,required=True)
    p.add_argument('--weight-mode',choices=('int_imf','power_law_initial_coordinate'),default='int_imf')
    p.add_argument('--alpha',type=float)
    p.add_argument('--escape-mass-tilt',type=float,default=0.,help='Sensitivity tilt of baseline weights in current mass; not an inferred escape mass function')
    p.add_argument('--age-years',type=float,required=True)
    p.add_argument('--metallicity-dex',type=float,required=True)
    p.add_argument('--source-description',required=True)
    p.add_argument('--synthetic-calibration',action='store_true')
    a=p.parse_args()
    pop=build_population(read_isochrone(a.isochrone),mass_basis=a.mass_basis,weight_mode=a.weight_mode,alpha=a.alpha,
                         age_years=a.age_years,metallicity_dex=a.metallicity_dex,
                         source_description=a.source_description,synthetic=a.synthetic_calibration,
                         escape_mass_tilt=a.escape_mass_tilt)
    pop.metadata.update(input_sha256=sha256_file(a.isochrone),input_path=str(a.isochrone.resolve()))
    pop.save(a.out)
    print(json.dumps(pop.metadata,indent=2))


if __name__=='__main__':
    main()
