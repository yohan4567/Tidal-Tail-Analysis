"""Copy and audit completed local/returned comparisons without changing live runs."""
import argparse
import json
import math
import re
import shutil
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp

from plot_mlr_parameters import HERE, STUDY, comparison_directories, document_tables, load_estimates, q_regions, read, sha

OFFSETS=(0.,-.5,-1.)
PARTS={0.:'offset_0',-.5:'offset_minus0p5',-1.:'offset_minus1'}
ALIASES={'NGC5904':'M5','NGC6121':'M4','NGC6171':'M107','NGC6341':'M92',
    'NGC7078':'M15','NGC7089':'M2','NGC7099':'M30','Pal5':'Palomar 5',
    'NGC6093':'M80','NGC6266':'M62','NGC6273':'M19','NGC5139':'ω Cen','NGC6656':'M22'}


def dump(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding='utf-8')


def copy_verified(source,dest,manifest,out):
    source,dest=Path(source),Path(dest)
    before=sha(source);dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists():raise FileExistsError('Do not overwrite a result snapshot: '+str(dest))
    shutil.copy2(source,dest)
    if before!=sha(dest) or before!=sha(source):raise ValueError('Source/copy changed during consolidation')
    manifest.append(dict(source=str(source.resolve()),destination=dest.relative_to(out).as_posix(),
                         sha256=before,bytes=dest.stat().st_size))


def verify_scores(path,offset):
    score=read(path);frame=pd.DataFrame(score['candidates'])
    if len(frame)!=1600 or frame.duplicated(['rate_msun_myr','duration_myr','scenario']).any():
        raise ValueError('Incomplete or duplicate raw score grid')
    if any(json.loads(v)['escape_slope_offset']!=offset for v in frame.scenario):raise ValueError('Mixed offsets')
    if not np.allclose(frame.q,2*(frame.negative_log_likelihood-frame.negative_log_likelihood.min()),atol=1e-9):
        raise ValueError('Raw score q mismatch')
    if path.with_suffix('.csv').exists():
        raw=pd.read_csv(path.with_suffix('.csv'))
        pd.testing.assert_frame_equal(raw[frame.columns],frame,check_dtype=False,rtol=1e-10,atol=1e-9)
    return score,frame


def run(out,doc):
    out=Path(out).resolve()
    if out.exists():raise FileExistsError('Choose a fresh consolidation snapshot')
    by_offset={d:load_estimates(d)[0] for d in OFFSETS}
    complete=set.intersection(*(set(v) for v in by_offset.values()))
    master_path=STUDY/'catalogs/master_clusters.csv';master=pd.read_csv(master_path).set_index('id')
    shapes,physical=document_tables(doc)
    if not complete<=set(physical):raise ValueError('Completed result outside the declared 29-cluster sample')
    for cid,p in physical.items():
        m=master.loc[cid]
        if not(abs(m['mass']-p['mass'])<1 and abs(m.ecc-p['ecc'])<=.0051 and abs(m.R_gc-p['R_gc'])<=.051):
            raise ValueError('Frozen physical catalogue and supplied table differ: '+cid)
    out.mkdir(parents=True);manifest=[];records=[];checks=[];grids=[];profiles=[];identities={}
    registry=read(HERE/'outputs/other_pc_all_remaining_20261101_v1/READY_PACKAGES.json')['ready_packages']
    old=read(HERE/'outputs/offset_results_20261101/summary.json')['results']
    for cid in sorted(complete,key=lambda c:(c=='Pal5',int(c[3:]) if c.startswith('NGC') else 99999)):
        for delta in OFFSETS:
            estimate=by_offset[delta][cid];part=PARTS[delta];dest=out/'completed'/cid/part
            dest.mkdir(parents=True)
            for p,h in estimate['sources'].items():
                if sha(p)!=h:raise ValueError('Source changed after selection: '+p)
            source_score=next(Path(p) for p in estimate['sources'] if p.endswith('joint_grid.json'))
            score,raw=verify_scores(source_score,delta)
            identities.setdefault(cid,set()).add(score.get('observation_identity'))
            returned=HERE.parent/'results' in source_score.parents
            folder=comparison_directories(delta).get(cid)
            if cid in ('NGC288','Pal5','NGC6171'):
                legacy=next(r for r in old if r['cluster']==cid and r['delta']==delta)
                grid=pd.DataFrame(legacy['likelihood_grid'])
                node_path=next(Path(p) for p in estimate['sources'] if p.endswith('integrated_nodes.csv'))
                nodes=pd.read_csv(node_path)
                prior_path=next(Path(p) for p in estimate['sources'] if p.endswith('population_prior.json'))
                integration_plan=read(next(Path(p) for p in estimate['sources'] if p.endswith('plan.json')))
                runtime_prior_verified=sha(prior_path)==integration_plan['population_prior_sha256']
                if not runtime_prior_verified:raise ValueError('Legacy population prior identity mismatch')
                for p in estimate['sources']:
                    p=Path(p)
                    if p==source_score:continue
                    sub='population' if p.name=='population_prior.json' else 'integration'
                    copy_verified(p,dest/'source'/sub/p.name,manifest,out)
            else:
                meta=read(folder/'comparison.json');grid=pd.read_csv(folder/'likelihood_grid.csv')
                node_path=folder/'integrated_nodes.jsonl';nodes=pd.read_json(node_path,lines=True)
                if returned:
                    prior_path=Path(registry[cid]['archive']).with_suffix('')/'inputs/populations'/part/'population_prior.json'
                else:prior_path=folder.parent/'population_bundle/population_prior.json'
                runtime_prior_verified=sha(prior_path)==meta['plan']['population_prior_sha256']
                if not returned and not runtime_prior_verified:raise ValueError('Local population prior identity mismatch')
                for p in sorted(folder.iterdir()):
                    if p.is_file():copy_verified(p,dest/'source/comparison'/p.name,manifest,out)
                copy_verified(prior_path,dest/'source/population/population_prior.json',manifest,out)
            prior=read(prior_path);weights=prior['scenario_weights']
            if not math.isclose(sum(weights.values()),1.,abs_tol=1e-12):raise ValueError('Population weights not normalized')
            if len(nodes)!=1600 or nodes.duplicated(['rate_msun_myr','duration_myr','scenario']).any():raise ValueError('Incomplete integrated nodes')
            keys=['rate_msun_myr','duration_myr','scenario']
            if set(map(tuple,nodes[keys].to_numpy()))!=set(map(tuple,raw[keys].to_numpy())):raise ValueError('Raw/integrated parameter nodes differ')
            calculated=[]
            for (rate,time),group in nodes.groupby(['rate_msun_myr','duration_myr']):
                if set(group.scenario)!=set(weights):raise ValueError('Population node weights do not match returned nodes')
                value=float(logsumexp([r.log_marginal_likelihood+np.log(weights[r.scenario]) for r in group.itertuples()]))
                calculated.append(dict(rate_msun_myr=float(rate),duration_myr=float(time),log_likelihood=value))
            expected=pd.DataFrame(calculated).sort_values(['rate_msun_myr','duration_myr']).reset_index(drop=True)
            expected['q']=2*(expected.log_likelihood.max()-expected.log_likelihood)
            grid=grid.sort_values(['rate_msun_myr','duration_myr']).reset_index(drop=True)
            pd.testing.assert_frame_equal(grid[expected.columns],expected,check_dtype=False,rtol=0,atol=1e-7)
            profile=expected.groupby('rate_msun_myr',as_index=False).q.min()
            supplied=pd.DataFrame(estimate['profile']).sort_values('rate_msun_myr').reset_index(drop=True)
            pd.testing.assert_frame_equal(profile,supplied,check_dtype=False,rtol=0,atol=1e-7)
            for p in sorted(source_score.parent.iterdir()):
                if p.is_file() and p.suffix in ('.json','.csv','.png','.md'):
                    copy_verified(p,dest/'source/scores'/p.name,manifest,out)
            expected.to_csv(dest/'likelihood_grid.csv',index=False,encoding='utf-8-sig')
            profile.to_csv(dest/'rate_profile.csv',index=False,encoding='utf-8-sig')
            expected.groupby('duration_myr',as_index=False).q.min().to_csv(dest/'time_profile.csv',index=False,encoding='utf-8-sig')
            best=estimate['best'];rate=float(best['rate_msun_myr']);duration=float(best['duration_myr'])
            x=profile.rate_msun_myr.to_numpy();q=profile.q.to_numpy();regions=q_regions(x,q)
            region=next(s for s in regions if s['lower']-1e-9<=rate<=s['upper']+1e-9)
            cat=master.loc[cid];mass=float(cat['mass'])
            row=dict(cluster=cid,name='NGC '+cid[3:] if cid.startswith('NGC') else 'Pal 5',alias=ALIASES.get(cid,''),
                delta=delta,origin='外機回傳' if returned else '本機',rate=rate,duration_myr=duration,
                q_lower=region['lower'],q_upper=region['upper'],lower_truncated=region['lower_grid_truncated'],
                upper_truncated=region['upper_grid_truncated'],diagnostic_minus=None if region['lower_grid_truncated'] else rate-region['lower'],
                diagnostic_plus=None if region['upper_grid_truncated'] else region['upper']-rate,
                q_regions=regions,rate_on_grid_boundary=rate in (float(x.min()),float(x.max())),time_on_grid_boundary=duration in (100.,2000.),
                mass_msun=mass,R_gc_kpc=float(cat.R_gc),eccentricity=float(cat.ecc),specific_rate_myr_inv=rate/mass,
                specific_rate_gyr_inv=1000*rate/mass,tail_shape=shapes[cid]['shape'],S1=shapes[cid]['S1'],S2=shapes[cid]['S2'],
                source_hashes=estimate['sources'],result_folder=str(dest),full_execution_provenance_verified=False)
            records.append(row);dump(dest/'result.json',row)
            g=expected.copy();g.insert(0,'cluster',cid);g.insert(1,'delta',delta);grids.append(g)
            rp=profile.copy();rp.insert(0,'cluster',cid);rp.insert(1,'delta',delta);profiles.append(rp)
            checks.append(dict(cluster=cid,delta=delta,source_score_sha256=sha(source_score),raw_grid_verified=True,
                copied_sources_verified=True,nuisance_integration_surface_recomputed=True,profile_recomputed=True,
                runtime_population_prior_hash_verified=runtime_prior_verified,
                population_weights_numeric_compatibility_verified=True,
                full_external_runtime_provenance_returned=False if returned else None,
                quadrature_max_relative_error=float(nodes.quadrature_relative_error.max())))
        if len(identities[cid])!=1:raise ValueError('Observation identity changes across offsets: '+cid)
    # Keep incomplete score-only returns outside the comparable result table.
    pending=[]
    for folder in sorted((HERE.parent/'results').iterdir()):
        if not folder.is_dir():continue
        cid='NGC'+folder.name.lower().replace('ngc','').strip()
        if cid in complete:continue
        for scorefile in sorted(folder.glob('*/scores/joint_grid.json')):
            score=read(scorefile);deltas={json.loads(r['scenario'])['escape_slope_offset'] for r in score['candidates']}
            if len(deltas)!=1:raise ValueError('Mixed offsets in pending return')
            delta=deltas.pop();score,raw=verify_scores(scorefile,delta)
            dest=out/'pending_comparison'/cid/PARTS[delta]/'scores'
            for p in sorted(scorefile.parent.iterdir()):
                if p.is_file():copy_verified(p,dest/p.name,manifest,out)
            pending.append(dict(cluster=cid,delta=delta,status='scores_only_not_mixed_into_integrated_results'))
    copy_verified(master_path,out/'catalogue/master_clusters.csv',manifest,out)
    now=datetime.now().astimezone().isoformat(timespec='seconds')
    base_status={}
    for line in (HERE/'outputs/cluster_status_20261101_latest.md').read_text(encoding='utf-8').splitlines():
        if line.startswith('| ') and not line.startswith('| 星團'):
            cells=[v.strip() for v in line.strip('|').split('|')]
            if len(cells)==3:
                m=re.search(r'NGC\s*(\d+)',cells[0]);cid='NGC'+m[1] if m else ('Pal5' if cells[0]=='Pal 5' else 'ESO452-SC11')
                base_status[cid]=(cells[1],cells[2])
    statuses=[]
    for cid in physical:
        cat=master.loc[cid];status,note=base_status.get(cid,('待確認','尚無完成的整合比較結果'))
        if cid in complete:
            origin=next(r['origin'] for r in records if r['cluster']==cid)
            status='已完成（'+origin+'）';note='三種 offset 的積分比較已核對並納入本次整合'
        statuses.append(dict(cluster=cid,alias=ALIASES.get(cid,''),status=status,note=note,
            included_in_completed_results=cid in complete,available_offsets=[d for d in OFFSETS if cid in by_offset[d]],
            mass_msun=float(cat['mass']),R_gc_kpc=float(cat.R_gc),eccentricity=float(cat.ecc),
            tail_shape=shapes[cid]['shape'],S1=shapes[cid]['S1'],S2=shapes[cid]['S2']))
    fields=['cluster','name','alias','origin','delta','rate','diagnostic_minus','diagnostic_plus','q_lower','q_upper',
        'lower_truncated','upper_truncated','rate_on_grid_boundary','duration_myr','time_on_grid_boundary',
        'mass_msun','R_gc_kpc','eccentricity','specific_rate_myr_inv','specific_rate_gyr_inv','tail_shape','S1','S2']
    frame=pd.DataFrame(records);frame[fields].to_csv(out/'mlr_all_offsets.csv',index=False,encoding='utf-8-sig')
    frame.loc[frame.delta==0,fields].to_csv(out/'mlr_offset0.csv',index=False,encoding='utf-8-sig')
    pd.concat(grids,ignore_index=True).to_csv(out/'likelihood_grids_all.csv',index=False,encoding='utf-8-sig')
    pd.concat(profiles,ignore_index=True).to_csv(out/'rate_profiles_all.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(statuses).drop(columns=['available_offsets']).to_csv(out/'cluster_status_29.csv',index=False,encoding='utf-8-sig')
    data=dict(checked_at=now,rows=records,statuses=statuses,pending=pending,
        cluster_count=len(complete),offset_result_count=len(records),no_target_posterior=True,
        confidence_intervals_calibrated=False,physical_rate_measurement_validated=False)
    dump(out/'data.json',data)
    dump(out/'source_manifest.json',dict(files=manifest,all_source_copies_sha256_verified=True))
    dump(out/'verification.json',dict(clusters=len(complete),offset_results=len(records),checks=checks,
        all_source_copies_sha256_verified=True,original_inputs_modified=False,confidence_intervals_calibrated=False,
        physical_rate_measurement_validated=False,score_only_returns_not_mixed=True,
        external_provenance_scope='Returned scores and comparison files checked; complete other-PC runtime and input provenance not returned. Local package weights reproduce the returned integrated surface, but cannot prove the unreturned runtime prior file identity.'))
    lines=['# 本機與外機結果整合','',f'資料快照：{now}。已整合 {len(complete)} 顆星團、{len(records)} 組 offset 結果。',
        '本資料夾是分析結果快照；原始 Gaia 全表、動力學粒子與執行環境仍保留在原專案，沒有移動或改寫。','',
        '## 檔案位置','',
        '- 星團結果總表.xlsx：MLR 總覽、診斷明細及全部 29 顆狀態。',
        '- mlr_all_offsets.csv、mlr_offset0.csv：完整精度數值，包含質量、銀心距離、偏心率及 MLR/M。',
        '- completed/：每顆星團、每種 offset 的統一比較網格、剖面、結果與原始來源副本。',
        '- pending_comparison/：僅有 scores 的外機回傳，另存而不混入積分比較。',
        '- source_manifest.json、verification.json：來源 SHA256、複製核對與數值檢查。','',
        '## 解讀','',
        'MLR 單位為 M☉/Myr，中央值是最佳已計算網格點。q≤1 診斷範圍在網格間線性插值，不是標準差或信賴區間；被網格截斷的方向留白。',
        'MLR/M 同時提供 Myr⁻¹ 與 Gyr⁻¹；只除以採用的現今目錄質量，未傳播分母的質量誤差。MLR 的人口質量範圍與目錄總質量不完全相同，因此它是條件式比較指標，不是已驗證的總質量分率。',
        'NGC6304、NGC6366 的三種 offset 最佳 MLR 均位於 31 的網格上界，不能當成已定位的精確流失率，也不能改寫成 >31 的物理下限。',
        '尾巴形態與 S1、S2 依使用者提供的資料表整理，未重新驗證形態偵測。質量與軌道數值採專案凍結目錄。',
        'NGC362 尚缺同一定義的 comparison，僅保留原始評分；其餘未完成或暫緩星團列在 29 顆狀態表中。',
        '舊 NGC288、Pal5、M107 的來源資料夾名稱含 posterior，但本次只重算背景與質量函數積分概似，未使用 R/T 先驗或後驗摘要。',
        '外機檔案已核對評分來源雜湊、完整網格、q、剖面及人口權重重建；尚未取得外機的全部執行來源，未宣稱完整科學驗證。']
    (out/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(dict(out=str(out),clusters=len(complete),offset_results=len(records),copied_files=len(manifest),
        returned_clusters=sorted(set(r['cluster'] for r in records if r['origin']=='外機回傳')),pending=pending),ensure_ascii=True))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--doc',type=Path,required=True);a=p.parse_args();run(a.out,a.doc)
