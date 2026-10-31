"""Maintain a sourced record for every target without inventing missing parameters."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from empirical_population import _Table,measured_constraint
from population_model import read_isochrone
from fetch_isochrone import validate_returned_model
from provenance import sha256_file
from run_cluster_survey import SAMPLE

HERE=Path(__file__).resolve().parent


def build(root):
    root=Path(root);adoptions=json.loads((HERE/'configs/survey_populations_v1.json').read_text())['clusters']
    def display_path(path):
        try:return str(Path(path).resolve().relative_to(HERE)).replace('\\','/')
        except ValueError:return str(Path(path).resolve())
    supplements=json.loads((HERE/'configs/population_supplements_v1.json').read_text())['clusters']
    master=pd.read_csv(root/'adopted_parameters.csv').set_index('id')
    html=HERE/'outputs/empirical_mf_source/parameter.html';records=[]
    for name,cid in SAMPLE:
        rule=dict(adoptions.get(cid,{}));rule.update(supplements.get(cid,{}))
        age=rule.get('age_gyr');row=master.loc[cid]
        finite=lambda v:float(v) if np.isfinite(v) else None
        try:mf=dict(status='measured_global_model_constraint',**measured_constraint(html,cid))
        except ValueError as exc:
            mf=dict(status='unavailable_for_direct_MF_quadrature',reason=str(exc))
            if cid=='ESO452-SC11':
                parser=_Table();parser.feed(html.read_text(encoding='utf-8'))
                rows=[r for r in parser.rows if len(r)>20 and re.search(r'ESO\s*452',r[0][0])]
                if len(rows)==1 and rows[0][20][1]:
                    nums=re.findall(r'[+-]?\d+(?:\.\d+)?',rows[0][20][0])
                    if len(nums)==2:mf.update(status='indirect_estimate_only',slope_mean=float(nums[0]),slope_sigma=float(nums[1]),
                        measured_mass_range_msun=None,source_sha256=sha256_file(html),source='https://people.smp.uq.edu.au/HolgerBaumgardt/globular/parameter.html')
        iso=(HERE/'outputs/Pal5_science_calibration/isochrone.dat') if cid=='Pal5' else root/cid/'calibration_inputs/isochrone.dat'
        imeta=iso.with_suffix('.json');iso_status='missing'
        if iso.exists() and imeta.exists():
            try:
                m=json.loads(imeta.read_text());assert m['sha256']==sha256_file(iso),'checksum mismatch'
                validate_returned_model(read_isochrone(iso),m['age_years'],m['metallicity_mh'])
                iso_status='validated_returned_age_and_metallicity'
            except (ValueError,KeyError,AssertionError) as exc:iso_status='invalid: '+str(exc)
        elif iso.exists():iso_status='unverified_model_without_metadata'
        blockers=[]
        if age is None:blockers.append('No approved single age; multi-population model required')
        if rule.get('hold_reason'):blockers.append(rule['hold_reason'])
        if mf['status']!='measured_global_model_constraint':blockers.append('No direct global MF with fitted mass interval')
        if not (root/cid/'selection/selection.json').exists():blockers.append('Missing real observed selection')
        if iso_status!='validated_returned_age_and_metallicity':blockers.append(iso_status)
        if finite(row.feh) is None:blockers.append('Frozen catalogue [Fe/H] missing; composition needs explicit adoption')
        if finite(row.ebv) is None:blockers.append('Frozen catalogue reddening missing')
        records.append(dict(cluster=cid,name=name,age=rule,global_mf=mf,
            composition=dict(frozen_catalogue_feh=finite(row.feh),frozen_catalogue_ebv=finite(row.ebv),
                source=display_path(root/'adopted_parameters.csv'),source_sha256=sha256_file(root/'adopted_parameters.csv'),
                mh_policy='Existing scaled-solar approximation; literature alpha/Fe recorded separately, never silently substituted'),
            isochrone=dict(path=display_path(iso),status=iso_status),
            escaping_population=dict(status='measured_tail_segment_with_extrapolation' if cid=='Pal5' else 'not_measured',
                delta_values=[0.,-.5,-1.],delta_status='sensitivity_assumptions_not_measurements',
                low_mass_extrapolation_constrained=False,remnant_escape_fraction=None,
                helium_distribution=None,binary_distribution=None,total_absolute_rate_supported=False),
            blockers=blockers))
    return dict(schema_version=1,updated_utc=datetime.now(timezone.utc).isoformat(),clusters=records,
        source_sha256=sha256_file(Path(__file__)),scope='Every target has a record; unknown values remain null, not invented')


def save(root):
    report=build(root)
    (HERE/'configs/population_registry_v1.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    lines=['# 29 顆星團人口參數清單','','有紀錄不等於參數已量測。δ、低質量外推與殘骸缺口仍須保留。','',
        '| 星團 | 年齡 Gyr | 全域 MF α | 等時線 | 缺口 |','|---|---|---|---|---|']
    for r in report['clusters']:
        age=r['age'].get('age_gyr');mf=r['global_mf'];alpha=mf.get('slope_mean')
        lines.append(f"| {r['name']} | {age if age is not None else '多人口；不可用單值替代'} | {alpha if alpha is not None else '缺'} ({mf['status']}) | {r['isochrone']['status']} | {'；'.join(r['blockers']) or '輸入初檢齊備；逃逸人口仍待約束'} |")
    lines+=['','所有參數來源、誤差、量測範圍及模型假設見 configs/population_registry_v1.json。',
        '年齡若為文獻模型採用值，不當成有誤差分布的直接量測；不為缺失誤差指定任意高斯。',
        'ω Cen 與已知金屬量複雜星團保留 hold，不為趕期限取消。']
    (HERE/'docs/POPULATION_REGISTRY_ZH.md').write_text('\n'.join(lines),encoding='utf-8')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=HERE/'outputs/cluster_survey_29_v1')
    a=p.parse_args();r=save(a.root);print('Recorded',len(r['clusters']),'clusters')
