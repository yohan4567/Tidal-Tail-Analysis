"""Measured global-MF constrained populations, retaining isochrone phase quadrature."""
from html.parser import HTMLParser
import json,re
from pathlib import Path
import numpy as np
from numpy.polynomial.hermite import hermgauss
from population_model import Population,build_population,read_isochrone
from provenance import sha256_file

class _Table(HTMLParser):
    def __init__(self):super().__init__();self.rows=[];self.row=None;self.cell=None;self.estimated=False
    def handle_starttag(self,tag,attrs):
        if tag=='tr':self.row=[]
        if tag in ('td','th'):self.cell=[];self.estimated=False
        if self.cell is not None:
            style=' '.join(v or '' for k,v in attrs if k in ('style','color')).lower()
            if 'green' in style or '#008000' in style or '#008888' in style:self.estimated=True
    def handle_data(self,data):
        if self.cell is not None:self.cell.append(data)
    def handle_endtag(self,tag):
        if tag in ('td','th') and self.cell is not None:
            if self.row is not None:self.row.append((' '.join(''.join(self.cell).split()),self.estimated))
            self.cell=None
        if tag=='tr' and self.row:self.rows.append(self.row);self.row=None


def measured_constraint(html,cluster):
    p=_Table();p.feed(Path(html).read_text(encoding='utf-8'))
    if cluster=='Pal5':pattern=r'\bPal\s*5\b'
    elif cluster=='ESO452-SC11':pattern=r'\bESO\s*452[\s-]*SC11\b'
    else:pattern=r'\bNGC\s*'+cluster.removeprefix('NGC')+r'\b'
    matches=[r for r in p.rows if len(r)>20 and re.search(pattern,r[0][0])]
    if len(matches)!=1:raise ValueError('Expected one official MF row for '+cluster)
    row=matches[0];text,estimated=row[20]
    if estimated:raise ValueError('Global MF is indirectly estimated, not a direct measured constraint: '+cluster)
    nums=re.findall(r'[+-]?\d+(?:\.\d+)?',text.replace('\u2212','-'))
    bounds=re.findall(r'\d+(?:\.\d+)?',row[19][0])
    if len(nums)!=2 or len(bounds)!=2:raise ValueError('Missing MF slope/error or fitted mass range')
    mu,sigma=map(float,nums);lo,hi=map(float,bounds)
    if not np.isfinite([mu,sigma,lo,hi]).all() or sigma<=0 or not 0<lo<hi:raise ValueError('Invalid MF constraint')
    return dict(cluster=cluster,slope_mean=mu,slope_sigma=sigma,convention='dN/dm proportional to m**slope',
        measured_mass_range_msun=[lo,hi],source='https://people.smp.uq.edu.au/HolgerBaumgardt/globular/parameter.html',
        paper='https://arxiv.org/abs/2303.01636',source_sha256=sha256_file(html),indirect_estimate=False,
        scope='Progenitor global present-day stellar MF, not a measured escaping MF')


def normal_nodes(mean,sigma,order=5):
    if not np.isfinite([mean,sigma]).all() or sigma<=0 or order<3:raise ValueError('Positive finite sigma and at least three quadrature nodes required')
    x,w=hermgauss(order)
    return mean+np.sqrt(2)*sigma*x,w/np.sqrt(np.pi)


def measured_population(isochrone,slope,*,age,mh,source):
    # Use the existing importer for phase/remnant exclusions and photometric nodes.
    base=build_population(read_isochrone(isochrone),mass_basis='isochrone_current_stellar_mass',
        source_description=source,age_years=age,metallicity_dex=mh)
    t=base.table.copy();m=t.initial_mass_msun.to_numpy();ms=t.loc[t.phase==1,['initial_mass_msun','current_mass_msun']]
    ms=ms.groupby('initial_mass_msun',sort=True).current_mass_msun.mean()
    if len(ms)<4 or (np.diff(ms.to_numpy())<0).any():raise ValueError('Need a monotone main-sequence initial/current mass map')
    xp=ms.index.to_numpy();yp=ms.to_numpy();proxy=np.interp(m,xp,yp)
    # For evolved phases, extrapolate the unevolved progenitor mapping. This
    # avoids treating a giant's wind-reduced current mass as a new birth mass.
    derivative=(yp[-1]-yp[-4])/(xp[-1]-xp[-4])
    if derivative<=0:raise ValueError('Invalid turnoff main-sequence mapping')
    high=m>xp[-1];proxy[high]=yp[-1]+derivative*(m[high]-xp[-1])
    if np.any(m<xp[0]-1e-7) or (np.diff(proxy)<0).any():raise ValueError('Unsupported pre-MS nodes or nonmonotone proxy')
    exponent=slope+1.;ratio=np.log(proxy[1:]/proxy[:-1])
    delta=ratio if abs(exponent)<1e-10 else proxy[:-1]**exponent*np.expm1(exponent*ratio)/exponent
    weights=.5*(np.r_[0.,delta]+np.r_[delta,0.])
    if (weights<0).any() or not weights.sum()>0:raise ValueError('Invalid MF quadrature')
    t['number_per_msun']=weights/np.sum(weights*t.current_mass_msun.to_numpy())
    if (t.number_per_msun<=0).any():raise ValueError('Nonpositive node support; cannot silently prune response nodes')
    meta=dict(base.metadata,weight_mode='empirical_present_day_ms_power_law_v1',escape_mass_tilt=None,
        population_builder_source_sha256=sha256_file(Path(__file__)),empirical_mf_slope=float(slope),source_description=source,input_sha256=sha256_file(isochrone),
        mass_function_status='Externally constrained progenitor MF plus an explicit escape-transfer scenario',
        escape_mass_function='dN/d(current main-sequence mass) proportional to mass**slope; evolved-phase counts retain initial-coordinate isochrone quadrature',
        represented_ms_mass_range_msun=[float(proxy.min()),float(proxy.max())],
        stars_per_basis_msun=float(t.number_per_msun.sum()),
        represented_current_mass_per_basis_msun=float(np.sum(t.current_mass_msun*t.number_per_msun)),
        rate_scope='Represented living-star mass only; excludes remnants/winds and stars outside supplied isochrone.',
        assumptions=['Single global-MF slope extrapolated outside its measured mass range.',
            'MS current/initial map is linearly extrapolated over narrow evolved progenitor interval; no independent evolved-phase IMF fit.',
            'Escape MF equality/offset relative to global GC MF is a scenario assumption, not an observation.',
            'Isochrone lower mass limit is retained explicitly; no invented photometry below model support.'])
    return Population(t,meta)


def prepare(cluster,isochrone,reference_config,html,out,*,offset=0.,order=5):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    constraint=measured_constraint(html,cluster);imeta=json.loads(Path(isochrone).with_suffix('.json').read_text())
    if imeta['sha256']!=sha256_file(isochrone):raise ValueError('Isochrone checksum changed')
    mean=constraint['slope_mean']+offset;nodes,weights=normal_nodes(mean,constraint['slope_sigma'],order)
    cfg=Path(reference_config).resolve();spec=json.loads(cfg.read_text());root=cfg.parent
    for k in ('quality','cache_dir','rate_history'):
        if spec.get(k):spec[k]=str((root/spec[k]).resolve())
    for k in ('path','metadata'):spec['errors'][k]=str((root/spec['errors'][k]).resolve())
    for k in ('path','metadata'):
        if spec['completeness'].get(k):spec['completeness'][k]=str((root/spec['completeness'][k]).resolve())
    configs=[];priors={}
    for i,(slope,weight) in enumerate(zip(nodes,weights)):
        source=f"{cluster}: measured global MF {constraint['slope_mean']} +/- {constraint['slope_sigma']}; escape slope offset {offset}; {constraint['paper']}; conditional extrapolation"
        pop=measured_population(isochrone,slope,age=imeta['age_years'],mh=imeta['metallicity_mh'],source=source)
        pop.metadata.update(global_mf_constraint=constraint,escape_slope_offset=offset,quadrature_prior_weight=float(weight))
        pop.save(out/f'population_{i}')
        scenario=dict(empirical_mf_slope=float(slope),escape_slope_offset=float(offset),parallax_zero_point_mas=spec['parallax_zero_point_mas'])
        s=dict(spec,population=str((out/f'population_{i}').resolve()),scenario=scenario,
            source_description=source+'; original cluster-specific quality/error/completeness calibration retained')
        s.pop('lookback_window_myr',None)
        path=out/f'response_{i}.json';path.write_text(json.dumps(s,indent=2),encoding='utf-8');configs.append(str(path.resolve()))
        priors[json.dumps(scenario,sort_keys=True)]=float(weight)
    report=dict(cluster=cluster,constraint=constraint,escape_slope_offset=offset,quadrature_order=order,
        scenario_weights=priors,response_configs=configs,source_sha256=sha256_file(Path(__file__)),
        reference_config_sha256=sha256_file(reference_config),isochrone_sha256=sha256_file(isochrone),
        interpretation='Gaussian quadrature of measured global-MF uncertainty; offset is a separate uncalibrated transfer scenario, not averaged with equal probabilities.')
    (out/'population_prior.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report
