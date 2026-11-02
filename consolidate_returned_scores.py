"""Audit and consolidate four externally returned score-only result sets."""
import argparse,hashlib,json,os
from pathlib import Path
import numpy as np
import pandas as pd

TARGETS=[('ngc362','NGC 362'),('6121','M4 / NGC 6121'),('6397','NGC 6397'),('7089','M2 / NGC 7089')]

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def run(source,out):
    source,out=Path(source).resolve(),Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    all_rows=[];best=[];checks=[];sources=[];profiles=[]
    for folder,cluster in TARGETS:
        found={}
        for path in (source/folder).glob('*/scores/joint_grid.csv'):
            data=pd.read_csv(path);meta=json.loads(path.with_suffix('.json').read_text())
            scenarios=data.scenario.map(json.loads)
            offsets={float(s['escape_slope_offset']) for s in scenarios}
            if len(offsets)!=1:raise ValueError('Mixed offsets in '+str(path))
            delta=next(iter(offsets))
            if delta in found:raise ValueError('Duplicate cluster/offset')
            found[delta]=path
            required=['duration_myr','rate_msun_myr','negative_log_likelihood','q','background_scale','expected_signal_stars']
            if not np.isfinite(data[required].to_numpy(float)).all():raise ValueError('Nonfinite scores')
            if data.duplicated(['duration_myr','rate_msun_myr','scenario']).any():raise ValueError('Duplicate parameter combinations')
            if len(data)!=1600 or data.rate_msun_myr.nunique()!=16 or data.duration_myr.nunique()!=20 or data.scenario.nunique()!=5:
                raise ValueError('Unexpected grid dimensions')
            if not all(len(part)==5 for _,part in data.groupby(['rate_msun_myr','duration_myr'])):raise ValueError('Incomplete nuisance nodes')
            nll=data.negative_log_likelihood.to_numpy();q=2*(nll-nll.min())
            if not np.allclose(q,data.q,rtol=1e-10,atol=1e-9):raise ValueError('q does not match the likelihood')
            original=pd.DataFrame(meta['candidates'])
            pd.testing.assert_frame_equal(data[original.columns],original,check_dtype=False,rtol=1e-10,atol=1e-9)
            for variable,name in [('rate_msun_myr','rate_profile.csv'),('duration_myr','time_profile.csv')]:
                supplied=pd.read_csv(path.parent/name).sort_values(variable).reset_index(drop=True)
                calculated=data.groupby(variable,as_index=False).q.min().sort_values(variable).reset_index(drop=True)
                pd.testing.assert_frame_equal(supplied,calculated,check_dtype=False,rtol=1e-10,atol=1e-9)
            row=data.loc[data.negative_log_likelihood.idxmin()];scenario=json.loads(row.scenario)
            minimum_rows=data.loc[np.isclose(data.negative_log_likelihood,row.negative_log_likelihood,rtol=0,atol=1e-10)]
            rate_lo,rate_hi=data.rate_msun_myr.min(),data.rate_msun_myr.max()
            time_lo,time_hi=data.duration_myr.min(),data.duration_myr.max()
            flags=[]
            if row.rate_msun_myr==rate_lo:flags.append('R 下界')
            if row.rate_msun_myr==rate_hi:flags.append('R 上界')
            if row.duration_myr==time_lo:flags.append('T 下界')
            if row.duration_myr==time_hi:flags.append('T 上界')
            bg=float(meta['background_only']['negative_log_likelihood'])
            rel=path.relative_to(source).as_posix();sha=digest(path)
            record=dict(cluster=cluster,delta=delta,rate_msun_myr=float(row.rate_msun_myr),duration_myr=float(row.duration_myr),
                lookback_gyr=float(row.duration_myr/1000),negative_log_likelihood=float(row.negative_log_likelihood),
                background_scale=float(row.background_scale),expected_signal_stars=float(row.expected_signal_stars),
                alpha_esc=float(scenario['empirical_mf_slope']),alpha_gc_node=float(scenario['empirical_mf_slope']-delta),
                background_only_nll=bg,two_delta_log_likelihood_vs_background=float(2*(bg-row.negative_log_likelihood)),
                grid_boundary='、'.join(flags) or '未落在邊界',best_tied_rows=len(minimum_rows),
                grid_points=320,population_nodes=5,score_rows=len(data),source_file=rel,source_sha256=sha,
                method='maximum_likelihood_profile_over_T_and_sampled_population',confidence_intervals_calibrated=False,
                physical_rate_measurement_validated=False,observation_identity=meta['observation_identity'])
            best.append(record)
            augmented=data.copy();augmented.insert(0,'cluster',cluster);augmented.insert(1,'delta',delta)
            augmented['alpha_esc']=[s['empirical_mf_slope'] for s in scenarios]
            augmented['alpha_gc_node']=augmented.alpha_esc-delta;augmented['source_file']=rel;all_rows.append(augmented)
            for variable in ['rate_msun_myr','duration_myr']:
                profile=data.groupby(variable,as_index=False).q.min();profile.insert(0,'cluster',cluster);profile.insert(1,'delta',delta)
                profile['profile_variable']=variable;profile=profile.rename(columns={variable:'value'});profiles.append(profile)
            checks.append(dict(cluster=cluster,delta=delta,score_rows=len(data),grid_complete=True,q_recomputed=True,
                json_csv_agree=True,supplied_profiles_agree=True,source_file=rel,sha256=sha))
            sources.append(dict(file=rel,sha256=sha,json_file=path.with_suffix('.json').relative_to(source).as_posix(),
                                json_sha256=digest(path.with_suffix('.json'))))
        if set(found)!={0.,-.5,-1.}:raise ValueError('Missing offset for '+cluster)
        identities={r['observation_identity'] for r in best if r['cluster']==cluster}
        if len(identities)!=1:raise ValueError('Observation identity differs across offsets')
    order={name:i for i,(_,name) in enumerate(TARGETS)}
    best.sort(key=lambda r:(order[r['cluster']],-r['delta']))
    df=pd.DataFrame(best);df.to_csv(out/'best_fits.csv',index=False,encoding='utf-8-sig')
    pd.concat(all_rows,ignore_index=True).to_csv(out/'all_model_scores.csv',index=False,encoding='utf-8-sig')
    pd.concat(profiles,ignore_index=True).to_csv(out/'profiles.csv',index=False,encoding='utf-8-sig')
    report=dict(clusters=4,offset_results=12,score_rows=sum(r['score_rows'] for r in checks),checks=checks,
        sources=sources,full_execution_provenance_verified=False,original_inputs_modified=False,
        compatible_with_nuisance_integrated_likelihood_results=False,
        missing_for_integrated_comparison='Original response/calibrated-count templates and population quadrature metadata; score-only CSVs do not retain background likelihood curvature.')
    (out/'validation.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    (out/'data.json').write_text(json.dumps(dict(best=best,validation=report),indent=2,ensure_ascii=False),encoding='utf-8')
    lines=['# 四組移機回傳結果統整','',
        '星團：NGC 362、M4（NGC 6121）、NGC 6397、M2（NGC 7089）。每顆三種 δ，共 12 組結果、19,200 筆模型評分。',
        'R 為所代表存活恆星的條件式平均逃逸質量率，單位 M☉/Myr。T 為回溯視窗，單位 Gyr。','',
        '| 星團 | R（δ=0） | R（δ=−0.5） | R（δ=−1） |','|---|---:|---:|---:|']
    for _,cluster in TARGETS:
        records=[r for r in best if r['cluster']==cluster]
        lines.append('| '+cluster+' | '+' | '.join(f"{r['rate_msun_myr']:g}" for r in records)+' |')
    lines+=['','## 定義與限制','',
        '本表直接使用回傳 scores：背景 β 在各點先最佳化，再對 T 與五個質量函數節點取最高概似。它與本機對背景、質量函數不確定性積分後的比較分數不同，不能混用成同一統計量。',
        '目前回傳檔沒有觀測模板與完整執行來源，無法僅由評分 CSV 重建上述積分。未計算 R/T 後驗，也未產生經校準的 68% 或 95% 信賴區間。',
        'M4 三種情境及 M2 的 δ=0 落在 R 下界 1。NGC 362 的 δ=−1 落在 R 上界 31；現有網格不能定位這些邊界外的最佳值。',
        'NGC 6397 與 M2 的最佳 T 均落在 0.1 Gyr 下界，不能視為精確回溯時間。T 不是星團年齡。',
        'M4 的最佳非零訊號模型在三種 δ 下都未優於同檔純背景模型，不能由這些最佳 R 推論已偵測到尾巴或非零流失。',
        '2Δℓ 相對純背景僅為分數改善，不是已校準的顯著性。q≤1 或 q≤3.84 不直接當成信賴區間。',
        '每顆星團才是一個研究樣本，三種 δ 是假設情境，不能當成三顆獨立星團。','',
        '已核對：每組 16×20×5 網格完整、參數組合未重複、CSV 與 JSON 相符、q 與原始概似相符、rate/time profile 可重算一致、同一星團三種 δ 的觀測識別碼一致。原始四組資料未修改。','',
        '## 各情境詳細數值','','| 星團 | δ | R（M☉/Myr） | T（Gyr） | 邊界 |','|---|---:|---:|---:|---|']
    for r in best:lines.append(f"| {r['cluster']} | {r['delta']:g} | {r['rate_msun_myr']:g} | {r['lookback_gyr']:g} | {r['grid_boundary']} |")
    (out/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(df[['cluster','delta','rate_msun_myr','lookback_gyr','grid_boundary']].to_string(index=False))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.source,a.out)
