"""Rescore frozen validation inputs with certified multimodal background profiling."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import joint_grid
import joint_recovery
from joint_fit import fixed_rate_fit
from provenance import sha256_file
from validate_joint_science import templates,mean_templates,compare


def finalize(root,out):
    root,out=Path(root).resolve(),Path(out).resolve()
    for stage in ('baseline','refined','independent_short','off_grid_truth'):
        if not (root/stage/'templates.json').exists(): raise ValueError('Raw dynamics/response stage must finish first: '+stage)
    if not (root/'injections.json').exists(): raise ValueError('Predeclared injection manifest missing')
    if out.exists(): raise ValueError('Use a fresh final analysis directory')
    # Raw response generation did not use this optimizer. Preserve raw source lineage;
    # only scores and bootstrap fits are replaced here, never particle data.
    joint_grid.fixed_rate_fit=fixed_rate_fit
    joint_recovery.fixed_rate_fit=fixed_rate_fit
    out.mkdir(parents=True)
    source_names=['finalize_joint_validation.py','joint_fit.py','joint_grid.py','joint_recovery.py']
    (out/'analysis_source_sha256.json').write_text(json.dumps({k:sha256_file(Path(__file__).parent/k) for k in source_names},indent=2))
    inputs={}
    for name in ('baseline','refined','independent_short','off_grid_truth'):
        paths=templates(root/name)
        for p in paths:
            meta=json.loads((p.parent/'response.json').read_text())
            if sha256_file(p)!=meta['template_sha256']: raise ValueError('Raw template integrity failure')
            inputs[str(p)]=sha256_file(p)
        joint_grid.analyse(paths,out/name,root/'regions.npz')
        print('Rescored with global optimizer:',name,flush=True)
    n,b,bmean,bgroups=mean_templates(templates(root/'baseline'))
    _,_,rmean,rgroups=mean_templates(templates(root/'refined'))
    _,_,smean,_=mean_templates(templates(root/'independent_short'))
    comparisons=[]
    for key,states in rmean.items():
        row=compare(bmean[key],states); row['duration_myr']=key[0]
        row.update(count_pass=row['relative_count_change']<=.05,shape_pass=row['spatial_total_variation']<=.1)
        comparisons.append(row)
    shortkey=next(iter(smean)); shortcheck=compare(rmean[shortkey],smean[shortkey])
    shortcheck.update(count_pass=shortcheck['relative_count_change']<=.05,shape_pass=shortcheck['spatial_total_variation']<=.1)
    recovery=joint_recovery.recover(templates(root/'baseline'),root/'injections.json',out/'recovery',100,seed=123)
    base=json.loads((out/'baseline/joint_grid.json').read_text()); ref=json.loads((out/'refined/joint_grid.json').read_text())
    passed=sum(d['noise'].get('major_pixel_target_passed',False) for d in base['diagnostics'])
    refpassed=sum(d['noise'].get('major_pixel_target_passed',False) for d in ref['diagnostics'])
    # Compare initial and corrected scoring to quantify impact on this dataset.
    old=pd.read_csv(root/'baseline/scores/joint_grid.csv'); new=pd.read_csv(out/'baseline/joint_grid.csv')
    merged=old.merge(new,on=['duration_myr','rate_msun_myr','scenario'],suffixes=('_old','_new'))
    delta=float(np.max(np.abs(merged.negative_log_likelihood_old-merged.negative_log_likelihood_new)))
    # Quantify age sensitivity without asserting a saturated regime.
    age_tables=[]
    for path in templates(root/'refined'):
        if '_T_1000' in path.parent.name:
            table=pd.read_csv(path.parent/'age_contributions.csv'); table['source']=str(path.parent); age_tables.append(table)
    age=pd.concat(age_tables,ignore_index=True)
    age.to_csv(out/'age_preflight.csv',index=False)
    fractions=[]
    for source,d in age.groupby('source'):
        total=d.expected_selected_stars.sum()
        fractions.append(dict(source=source,older_than_750_myr_fraction=float(d.loc[d.lookback_start_myr>=750,'expected_selected_stars'].sum()/total) if total else None))
    checks=dict(stratified_mass_normalization=True,age_response_accounting=True)
    for folder in ('baseline','refined','independent_short','off_grid_truth'):
        for path in (root/folder).glob('rate_*/simulation/particles.csv'):
            frame=pd.read_csv(path); config=json.loads((path.parent/'run.json').read_text())['config']
            checks['stratified_mass_normalization'] &= bool(np.isclose(frame.mass_weight_msun.sum(),config['duration_myr']*config['reference_mdot_msun_myr']))
        for path in templates(root/folder):
            meta=json.loads((path.parent/'response.json').read_text())
            ages=pd.read_csv(path.parent/'age_contributions.csv')
            checks['age_response_accounting'] &= bool(np.isclose(ages.expected_selected_stars.sum(),meta['expected_stars_in_fit_region']))
    result=dict(completed=True,physical_rate_measurement_validated=False,checks=checks,
        baseline_mc_target_passed=passed,baseline_candidates=len(base['diagnostics']),
        refined_mc_target_passed=refpassed,refined_candidates=len(ref['diagnostics']),
        density_comparisons=comparisons,independent_short_comparison=shortcheck,
        corrected_optimizer_maximum_nll_change=delta,age_sensitivity=fractions,
        recovery=recovery,input_sha256=inputs,
        raw_source_archive=str(root/'source_used'),
        caveats=['A completed validation experiment can fail its convergence targets.',
            'No independent N-body, full selection refitting or total population-mass validation.',
            'No observed physical MLR confidence interval is certified.'])
    (out/'validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    pd.DataFrame(comparisons).to_csv(out/'refinement.csv',index=False)
    lines=['# Pal5 聯合推論實測與驗證結果','',
        '這份報告完成了條件式數值與統計實驗；完成執行不等於所有檢查通過，也不等於取得已驗證的物理流失率。','',
        '## 實驗範圍','',
        '實際 Gaia Pal5 選星與校正輸入；MLR 4、9、25 M☉/Myr；T 250、500、750、1000 Myr；三個人口情境；每組三個獨立動力學種子。基準每個最長時間模型 400 粒子，加密對照為 800 粒子。另有獨立短時間模型與 13 M☉/Myr 的非網格真值。','',
        '## 數值結果','',
        f"基準逐像素 MC 目標通過 {passed}/{len(base['diagnostics'])}；加密模型通過 {refpassed}/{len(ref['diagnostics'])}。這是保守目標，未通過者不能宣稱模板已收斂。",'',
        '| T Myr | 加密後計數變化 | 空間 TV | 計數目標通過 | 形態目標通過 |',
        '|---|---|---|---|---|']
    for row in comparisons:
        lines.append(f"| {row['duration_myr']:g} | {row['relative_count_change']:.2%} | {row['spatial_total_variation']:.3f} | {row['count_pass']} | {row['shape_pass']} |")
    lines+=['','獨立短時間比較：', '',json.dumps(shortcheck,indent=2),'',
        '## 統計修正與回收測試','',
        '修正多校正狀態的背景 profile 搜尋，使用每個 Poisson 分量的凹性建立全域界限，避免單峰搜尋漏解。舊來源及舊分數保留，這裡重新計算分數與回收測試。',
        f'本資料修正前後最大 NLL 差異為 {delta:.6g}；極端案例回歸測試仍有必要。','',
        '六個真值情境各 100 次，前半校準、後半獨立檢查。涵蓋率及假陽性附 Wilson 區間。非網格情境不宣稱已校準連續區間。','',
        '| 情境 | rate bias | held-out coverage | held-out false positive |',
        '|---|---|---|---|']
    for case in recovery['scenarios']:
        lines.append(f"| {case['name']} | {case.get('rate_bias','—')} | {case.get('coverage','—')} | {case.get('held_out_false_positive_fraction','—')} |")
    lines+=['',f"尚未驗證的情境標記：{', '.join(recovery['missing_required_scenario_tags']) or '無'}。時間標記需要年代貢獻支持，不能為了補齊標記而假設已飽和。",'',
        '## 使用結論','',
        '程式流程、權重及條件式推論已實測；目前未將模板收斂與人口等系統誤差全部驗證通過，因此最佳 MLR 仍只作條件式試驗值。後續需增加取樣或在經驗證的空間尺度上比較，並檢查人口、背景與銀河勢。','',
        '[聯合圖](baseline/joint_grid.png) · [完整數值](validation.json) · [回收細節](recovery/recovery.json)']
    (out/'REPORT_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('Final corrected validation completed:',out,flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--raw',type=Path,required=True); p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); finalize(a.raw,a.out)
