"""Versioned literature/population repair of M92, M30, M15 and ESO452-SC11."""
import argparse,importlib.util,json,os,re,shutil,subprocess,sys,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from provenance import sha256_file
from transfer_status_io import write_json
from fetch_isochrone import fetch,validate_returned_model
from population_model import read_isochrone

HERE=Path(__file__).resolve().parent
POLICY=HERE/'configs/four_population_repairs_v1.json'
REVIEW=HERE/'outputs/other_pc_review_nine_20261101_v1/review_work_package'
OLD=HERE/'outputs/input_repairs_20261101_v1'

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,d):write_json(p,d,status_file=Path(p).name=='status.json')

def global_metallicity(feh,alpha,a=.638):
    if not np.isfinite([feh,alpha,a]).all() or not 0<a<1:raise ValueError('Finite abundances and a metal fraction in (0,1) required')
    return float(feh+np.log10(a*10**alpha+1-a))

def builder():
    code=OLD/'frozen_code';manifest=read(OLD/'repair_manifest.json')
    for name in ('empirical_population.py','population_model.py','population_quadrature_repair.py'):
        p=(code/name).resolve()
        if sha256_file(p)!=manifest['frozen_sha256'][str(p)]:raise ValueError('Reviewed population builder changed')
        if name!='empirical_population.py' and sha256_file(HERE/name)!=sha256_file(p):raise ValueError('Population dependency differs from reviewed snapshot')
    spec=importlib.util.spec_from_file_location('_four_repair_fixed_population',code/'empirical_population.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def obtain(case,age,mh,source,cached=None):
    case.mkdir(parents=True,exist_ok=True)
    for meta in sorted(case.glob('isochrone_attempt*.json')):
        data=read(meta);p=meta.with_suffix('.dat')
        if p.exists() and sha256_file(p)==data['sha256']:
            validate_returned_model(read_isochrone(p),age,mh);return p
    i=len(list(case.glob('isochrone_attempt*.dat')))+1;p=case/f'isochrone_attempt{i}.dat'
    if cached:
        data=read(cached.with_suffix('.json'))
        if sha256_file(cached)!=data['sha256']:raise ValueError('Cached isochrone checksum mismatch')
        validate_returned_model(read_isochrone(cached),age,mh)
        shutil.copy2(cached,p);shutil.copy2(cached.with_suffix('.json'),p.with_suffix('.json'))
    else:fetch(age,mh,p,source=source)
    return p

def cmd_review(iso,selection,ebv,out):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    table=read_isochrone(iso);r=read(selection/'selection.json');ridge=r['cmd_ridge']
    distance=float(r['cluster_selection_parameters']['R_sun']);dm=5*np.log10(distance)+10
    ag,abp,arp=np.array([.83627,1.08337,.63439])*3.1*ebv
    g=table.Gmag.to_numpy()+dm+ag;colour=(table.G_BPmag-table.G_RPmag).to_numpy()+abp-arp
    rg=np.asarray(ridge['G']);rc=np.asarray(ridge['color']);width=np.asarray(ridge['width']);delta=np.full(len(rg),np.inf)
    for phase in (1,2,3):
        mask=(table.label.to_numpy()==phase)&np.isfinite(g)&np.isfinite(colour)
        if mask.sum()<2:continue
        order=np.argsort(g[mask]);gg,cc=g[mask][order],colour[mask][order]
        gg,unique=np.unique(gg,return_index=True);cc=cc[unique]
        inside=(rg>=gg.min())&(rg<=gg.max())
        delta[inside]=np.minimum(delta[inside],np.abs(np.interp(rg[inside],gg,cc)-rc[inside]))
    supported=np.isfinite(delta);n=int(supported.sum())
    report=dict(cluster=r['cluster'],isochrone_sha256=sha256_file(iso),selection_sha256=sha256_file(selection/'selection.json'),
        reference_distance_kpc=distance,ebv=ebv,extinction_approximation='Constant G2V coefficients, Rv=3.1; temperature-dependent/differential extinction not modelled.',
        ridge_points=len(rg),supported_ridge_points=n,
        median_absolute_bp_rp_offset=float(np.median(delta[supported])) if n else None,
        fraction_within_empirical_ridge_width=float(np.mean(delta[supported]<=width[supported])) if n else None,
        scope='Diagnostic at cluster distance for MS/SGB/RGB; not a calibrated goodness-of-fit test or approval of the full stellar population. HB, dust, completeness and mass-function assumptions need separate checks.')
    dump(out.with_suffix('.json'),report)
    fig,ax=plt.subplots(figsize=(5,6));keep=(table.label.to_numpy()!=9)&(g<22)
    ax.scatter(colour[keep],g[keep],s=7,label='Adopted isochrone',alpha=.7)
    ax.errorbar(rc,rg,xerr=width,fmt='o',ms=3,color='black',label='Frozen empirical CMD ridge')
    ax.set(xlabel='Gaia BP - RP (mag)',ylabel='Gaia G (mag)',title=r['cluster']);ax.invert_yaxis();ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out.with_suffix('.png'),dpi=160);plt.close(fig)
    return report

def stage_study(root,clusters):
    study=root/'study';(study/'catalogs').mkdir(parents=True,exist_ok=True);(study/'gaia_fields').mkdir(exist_ok=True)
    for name in ('step3_tail_search.py','skymath.py'):
        if not (study/name).exists():shutil.copy2(REVIEW/'study'/name,study/name)
    catalog=study/'catalogs/master_clusters.csv'
    if not catalog.exists():
        frame=pd.read_csv(REVIEW/'study/catalogs/master_clusters.csv')
        frame.loc[frame.id=='ESO452-SC11','ebv']=.54
        frame.loc[frame.id=='ESO452-SC11','feh']=-.88
        frame.to_csv(catalog,index=False)
        dump(study/'catalogue_changes.json',dict(source_sha256=sha256_file(REVIEW/'study/catalogs/master_clusters.csv'),
            changed_cluster='ESO452-SC11',changes=dict(ebv=.54,feh=-.88),source='https://arxiv.org/abs/2502.01741v2',
            policy='New study copy only. Original membership/orbit distance and frozen selection are retained; [M/H]=-0.80 remains a separate population input.'))
    for cid in clusters:
        original=REVIEW/'inputs'/cid/'observed';dest=root/'observed'/cid
        if not dest.exists():
            dest.mkdir(parents=True);shutil.copytree(original/'selection',dest/'selection')
            shutil.copy2(original/'adopted_dynamics.json',dest/'adopted_dynamics.json')
        s=read(dest/'selection/selection.json')
        if s['selection_file_sha256']!=sha256_file(study/'step3_tail_search.py'):raise ValueError('Selector changed')
        field=REVIEW/'study/gaia_fields'/(cid+'_field.csv.gz')
        if sha256_file(field)!=s['source_field_sha256']:raise ValueError('Observed field changed')
        for suffix in ('_field.csv.gz','_field.json','_missing_tiles.json'):
            p=field.parent/(cid+suffix)
            if p.exists() and not (study/'gaia_fields'/p.name).exists():shutil.copy2(p,study/'gaia_fields'/p.name)
    return study

def stage_code(root,study):
    base=root/'model_base';code=base/'code'
    if (base/'repair_manifest.json').exists():return base,code
    original=read(OLD/'repair_manifest.json')
    for name,sha in original['frozen_sha256'].items():
        if sha256_file(name)!=sha:raise ValueError('Source snapshot changed')
    shutil.copytree(original['snapshot'],code,ignore=shutil.ignore_patterns('__pycache__','outputs'))
    for name in ('bright_calibration.py','quality_repair.py','likelihood_comparison.py','reviewed_population_policy.py'):
        shutil.copy2(HERE/name,code/name)
    p=code/'bright_calibration.py';text=p.read_text(encoding='utf-8')
    start=text.index("    if min(quality.arrays['n_parent'].min(),quality.arrays['n_available'].min())<quality.minimum:")
    end=text.index('    # The original MF table',start);text=text[:start]+text[end:]
    text=text.replace('from empirical_population import prepare as population','from reviewed_population_policy import prepare as population')
    text=text.replace('offset=offset,order=5)','offset=offset,order=5,allow_indirect_model=(cluster==\'ESO452-SC11\'))')
    text=text.replace("compression='gzip'","compression={'method':'gzip','compresslevel':1}");p.write_text(text,encoding='utf-8')
    p=code/'prepare_science_calibration.py';text=p.read_text(encoding='utf-8')
    text=text.replace('allsky_uniq_10.fits','allsky_M10_hpx7.hdf5').replace("mode='multi'","mode='hpx7'")
    text=text.replace("compression='gzip'","compression={'method':'gzip','compresslevel':1}");p.write_text(text,encoding='utf-8')
    shutil.copy2(OLD/'mf_source.html',base/'mf_source.html')
    files=[*code.glob('*.py'),*list((code/'configs').rglob('*.json')),base/'mf_source.html',POLICY]
    dump(base/'repair_manifest.json',dict(source=str(root/'observed'),snapshot=str(code),study=str(study),
        scope='Four-cluster conditional composition repair; ESO indirect MF requires explicit opt-in.',
        preparation_source_sha256=sha256_file(Path(__file__)),frozen_sha256={str(p.resolve()):sha256_file(p) for p in files}))
    return base,code

def verified_parent(cid,selection,target):
    from selection_model import parent_query,download_parent
    from adaptive_parent import fetch as adaptive_fetch
    s=read(selection/'selection.json');ra,dec=s['center_ra_deg'],s['center_dec_deg'];radius=s['field_radius_deg']+.05
    if target.exists():
        m=read(target.with_suffix('.json'))
        if m.get('query_status')!='OK' or m['sha256']!=sha256_file(target):raise ValueError('Invalid parent checkpoint')
        return
    if target.name=='parent_random10.csv.gz':
        cached=HERE/'outputs/nine_cluster_real_inputs_v1'/cid/target.name
        if cached.exists() and cached.with_suffix('.json').exists():
            m=read(cached.with_suffix('.json'))
            if m['sha256']==sha256_file(cached) and m.get('query_status')=='OK' and m.get('g_min')==12 and m.get('g_max')==21 and m.get('random_sampling_modulus')==10 and m.get('radius_deg')==radius and np.allclose(m.get('center'),[ra,dec]):
                shutil.copy2(cached,target);shutil.copy2(cached.with_suffix('.json'),target.with_suffix('.json'));return
    bright=target.name=='bright_parent.csv.gz';gmin,gmax,mod=(0.,12.,1) if bright else (12.,21.,10)
    query=parent_query(ra,dec,radius,gmin,gmax,mod);last=None
    for endpoint in ('https://gaia.aip.de/tap/sync','https://gaia.aip.de/tap/sync','https://gea.esac.esa.int/tap-server/tap/sync'):
        try:
            m=download_parent(query,target,100000 if bright else 300000,endpoint)
            m.update(center=[ra,dec],radius_deg=radius,g_min=gmin,g_max=gmax,random_sampling_modulus=mod)
            dump(target.with_suffix('.json'),m);return
        except RuntimeError as exc:
            if 'OVERFLOW' not in str(exc) or bright:raise
            adaptive_fetch(target,ra,dec,radius,endpoint='https://gaia.aip.de/tap/sync');return
        except OSError as exc:
            last=exc
            if target.exists() or target.with_suffix('.json').exists():raise
            time.sleep(1)
    raise RuntimeError('Parent endpoints unavailable; preserve query history and retry separately') from last

def run(root,inputs_only=False,selected=None):
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True);policy=read(POLICY)
    chosen=selected or list(policy['clusters'])
    if not chosen or len(set(chosen))!=len(chosen) or any(cid not in policy['clusters'] for cid in chosen):
        raise ValueError('Choose distinct reviewed clusters from the composition policy')
    clusters={cid:policy['clusters'][cid] for cid in chosen}
    states=read(root/'status.json').get('clusters',{}) if (root/'status.json').exists() else {}
    def update(cid,stage,**extra):
        states.setdefault(cid,{}).update(stage=stage,updated_utc=datetime.now(timezone.utc).isoformat(),physical_rate_measurement_validated=False,**extra)
        dump(root/'status.json',dict(clusters=states));print(cid,stage,extra.get('error',''),flush=True)
    for cid in clusters:
        if states.get(cid,{}).get('stage')!='conditional_model_package_ready':update(cid,'queued')
    fixed=builder();study=stage_study(root,clusters);os.environ['GC_TIDAL_STUDY']=str(study)
    base,code=stage_code(root,study)
    def headroom(cid):
        from run_deadline_queue import memory_snapshot
        while True:
            m=memory_snapshot()
            if not m or m['available_gib']>=5:return
            update(cid,'waiting_for_RAM',available_gib=m['available_gib']);time.sleep(30)
    def command(script,args,cid):
        headroom(cid)
        with (root/(cid+'_model.log')).open('a',encoding='utf-8') as log:
            p=subprocess.run([sys.executable,'-u',str(script),*map(str,args)],cwd=HERE,stdout=log,stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
        if p.returncode:raise RuntimeError('\n'.join((root/(cid+'_model.log')).read_text(encoding='utf-8',errors='replace').splitlines()[-10:]))
    catalog=pd.read_csv(study/'catalogs/master_clusters.csv').set_index('id')
    for cid,adoption in clusters.items():
        if states.get(cid,{}).get('stage')=='conditional_model_package_ready':continue
        try:
            dest=root/'observed'/cid;inp=dest/'calibration_inputs';inp.mkdir(exist_ok=True)
            cases=[('nominal',adoption.get('adopted_alpha_fe'))]
            if 'sensitivity_alpha_fe' in adoption:cases.append(('alpha_0p30',adoption['sensitivity_alpha_fe']))
            reviews=[]
            for label,alpha in cases:
                mh=adoption['metallicity_mh'] if alpha is None else global_metallicity(adoption['population_feh'],alpha)
                update(cid,'obtaining_verified_isochrone',case=label,target_mh=mh)
                cached=HERE/'outputs/cluster_survey_29_v1'/cid/'calibration_inputs/isochrone.dat' if alpha is None else None
                iso=obtain(root/'cases'/cid/label,adoption['age_gyr']*1e9,mh,
                    f'{cid}; {adoption["source"]}; literature-adopted composition; alpha={alpha}; classical Salaris surrogate, not direct alpha-enhanced tracks; case={label}',cached)
                params=dict(adoption,case=label,target_mh=mh,adopted_alpha_fe=alpha,
                    conversion=policy['composition_conversion'] if alpha is not None else None,
                    original_selection_distance_kpc=read(dest/'selection/selection.json')['cluster_selection_parameters']['R_sun'],
                    scientific_scope=policy['interpretation'],sensitivity_policy=policy['sensitivity_policy'])
                dump(iso.parent/'population_adoption.json',params)
                ebv=float(catalog.loc[cid,'ebv']);review=cmd_review(iso,dest/'selection',ebv,iso.parent/'cmd_review');reviews.append(dict(case=label,**review))
                if label=='nominal':
                    for suffix in ('.dat','.json'):shutil.copy2(iso.with_suffix(suffix),inp/('isochrone'+suffix))
                    for name in ('population_adoption.json','cmd_review.json','cmd_review.png'):shutil.copy2(iso.parent/name,inp/name)
                else:
                    for suffix in ('.dat','.json'):shutil.copy2(iso.with_suffix(suffix),inp/('isochrone_'+label+suffix))
                constraint=fixed.measured_constraint(base/'mf_source.html',cid) if cid!='ESO452-SC11' else __import__('reviewed_population_policy').indirect_constraint(base/'mf_source.html',cid)
                for offset in (0.,-.5,-1.):
                    nodes,weights=fixed.normal_nodes(constraint['slope_mean']+offset,constraint['slope_sigma'],5)
                    for i,(slope,weight) in enumerate(zip(nodes,weights)):
                        pop=fixed.measured_population(iso,slope,age=adoption['age_gyr']*1e9,mh=mh,source=f'{cid}: reviewed {label} population; MF kind '+('indirect model assumption' if cid=='ESO452-SC11' else 'measured global MF'))
                        if not (pop.table.number_per_msun>0).all() or not np.isclose(np.sum(pop.table.number_per_msun*pop.table.current_mass_msun),1.,rtol=1e-12):raise ValueError('Population normalization/support failed')
                        pop.metadata.update(global_mf_constraint=constraint,escape_slope_offset=offset,quadrature_prior_weight=float(weight),composition_adoption=params)
                        folder=iso.parent/'population_checks'/f'offset_{offset:g}'/f'population_{i}'
                        if not folder.exists():pop.save(folder)
            dump(root/'cases'/cid/'review.json',dict(cmd_diagnostics=reviews,population_nodes_per_case=15,absolute_rate_validated=False))
            cfg=read(dest/'adopted_dynamics.json');cfg['mass_evolution']=dict(mode='linear_escape',dynamical_to_population_mass=1.,population_mass_basis='isochrone_current_stellar_mass',source_description='Conditional represented living-star baseline f=1; fixed Plummer scale; no remnant/wind correction.')
            dump(dest/'coupled_dynamics.json',cfg)
            update(cid,'isochrone_and_population_inputs_checked',cases=len(cases),mass_function_indirect=cid=='ESO452-SC11')
            if inputs_only:continue
            headroom(cid);update(cid,'preparing_verified_gaia_parents')
            verified_parent(cid,dest/'selection',inp/'parent_random10.csv.gz');verified_parent(cid,dest/'selection',inp/'bright_parent.csv.gz')
            extension=root/'model_repairs'/cid
            if not (extension/'bright_repair.json').exists():command(code/'bright_calibration.py',['--repair',base,'--cluster',cid,'--bright',inp/'bright_parent.csv.gz','--out',extension],cid)
            for name in ('population_adoption.json','cmd_review.json','cmd_review.png'):
                shutil.copy2(inp/name,extension/cid/'calibration_inputs'/name)
            for p in inp.glob('isochrone_alpha*'):shutil.copy2(p,extension/cid/'calibration_inputs'/p.name)
            if not (extension/'preflight.json').exists():
                update(cid,'testing_real_input_response');command(code/'bright_calibration.py',['--finalize','--out',extension],cid)
            package=root/'packages'/(cid+'_portable_v1')
            if not package.with_suffix('.zip').exists():
                update(cid,'exporting_conditional_package');headroom(cid)
                from portable_cluster import export
                export(extension,cid,package,analysis_kind='model_comparison')
            validation=root/'package_validation'/cid
            if not (validation/'validation.json').exists():command(HERE/'validate_transfer_package.py',['--bundle',package,'--simulation',extension/'preflight/simulation','--out',validation],cid)
            update(cid,'conditional_model_package_ready',archive=str(package.with_suffix('.zip')),mass_function_indirect=cid=='ESO452-SC11',
                note='Verified software/population inputs and relocated small response. CMD diagnostics and population assumptions remain conditional; no validated absolute-rate claim.')
        except Exception as exc:update(cid,'repair_or_calibration_failed',error=str(exc))
    dump(root/'finished.json',dict(clusters=states,finished_utc=datetime.now(timezone.utc).isoformat(),physical_rate_measurement_validated=False))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--inputs-only',action='store_true')
    p.add_argument('--clusters',nargs='+',choices=('NGC6341','NGC7099','NGC7078','ESO452-SC11'));a=p.parse_args()
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
    run(a.root,a.inputs_only,a.clusters)
