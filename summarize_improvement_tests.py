"""Summarize dense trajectories, age preflight and independent profile calibration."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from profile_recovery import rank_threshold
from joint_recovery import wilson
from provenance import sha256_file

HERE=Path(__file__).resolve().parent


def run():
    root=HERE/'outputs';out=root/'Pal5_improvement_tests_v1';out.mkdir(exist_ok=True)
    dense=root/'Pal5_dense_8000_analysis_v1'
    noise=pd.read_csv(dense/'noise.csv');nested=pd.read_csv(dense/'nested_thinning.csv')
    pairs=pd.read_csv(dense/'independent_seed_pairs.csv')
    profile=pd.read_csv(dense/'conditional_time_profiles.csv')
    ages=pd.read_csv(root/'Pal5_age_preflight_4gyr_v1/age_count_preflight.csv')
    old=pd.read_csv(root/'Pal5_profile_recovery_v2/trials.csv')
    old=old[old.name=='exact_mean_T500']
    new=pd.read_csv(root/'Pal5_profile_calibration_independent/trials.csv')
    newmeta=json.loads((root/'Pal5_profile_calibration_independent/profile_recovery.json').read_text())
    if len(old)!=500 or len(new)!=1000 or newmeta['seed']!=271828:raise ValueError('Wrong follow-up bank or seed')
    enlarged=[]
    for target in ('rate_profile','time_profile','time_rate_joint','candidate_point'):
        calibration=pd.concat([old[target],new.loc[new.partition=='calibration',target]])
        held=new.loc[new.partition=='held_out',target]
        threshold,k=rank_threshold(calibration);covered=int((held<=threshold).sum())
        enlarged.append(dict(target=target,calibration_n=len(calibration),rank=k,validation_n=len(held),
            threshold=threshold,covered=covered,coverage=covered/len(held),wilson_95=wilson(covered,len(held))))
    original=json.loads((root/'Pal5_profile_recovery_v2/profile_recovery.json').read_text())
    fraction=(ages.loc[ages.T_myr==4000,'mean'].iloc[0]-ages.loc[ages.T_myr==1000,'mean'].iloc[0])/ages.loc[ages.T_myr==4000,'mean'].iloc[0]
    shapes=[]
    for factor in (1,2,4,8,16):
        p=profile[(profile.factor==factor)&(profile.seed=='mean')].pivot(index='T_myr',columns='density',values='q')
        shapes.append(dict(factor=factor,max_nested_q_change_4k_8k=float(abs(p[4000]-p[8000]).max())))
    files=[dense/'convergence.json',dense/'noise.csv',dense/'nested_thinning.csv',dense/'independent_seed_pairs.csv',
        dense/'conditional_time_profiles.csv',root/'Pal5_age_preflight_4gyr_v1/age_count_preflight.csv',
        root/'Pal5_profile_recovery_v2/profile_recovery.json',root/'Pal5_profile_recovery_v2/trials.csv',
        root/'Pal5_profile_calibration_independent/profile_recovery.json',root/'Pal5_profile_calibration_independent/trials.csv',
        out/'unittest.log']
    testlog=(out/'unittest.log').read_text(encoding='utf-8',errors='replace')
    if 'Ran 60 tests' not in testlog or '\nOK' not in testlog:raise ValueError('Expected successful full test log')
    result=dict(completed=True,unit_tests_passed=60,dense_tracers=24000,age_preflight_tracers=1200,
        independent_scored_datasets=2500,physical_rate_measurement_validated=False,
        dense_scope='Rate9/population0, three seeds; no high-density full MLR grid.',
        fine_pixels_mc_passed=int(noise.loc[noise.factor==1,'major_pixel_target_passed'].sum()),
        coarse_factor4_mc_passed=int(noise.loc[noise.factor==4,'major_pixel_target_passed'].sum()),
        fine_independent_pair_tv_range=[float(pairs.loc[pairs.factor==1,'TV'].min()),float(pairs.loc[pairs.factor==1,'TV'].max())],
        fine_nested_4k_8k_tv_range=[float(nested.loc[(nested.factor==1)&(nested.lower_density==4000),'TV'].min()),float(nested.loc[(nested.factor==1)&(nested.lower_density==4000),'TV'].max())],
        profile_stability=shapes,age_preflight_after_1gyr_fraction=float(fraction),
        enlarged_calibration=enlarged,original_profile_recovery=original,
        input_sha256={str(p.resolve()):sha256_file(p) for p in files},source_sha256=sha256_file(__file__),
        limitations=['MC variances estimated from only three dynamics seeds.',
            'Same observation sampling stream/quality states across dynamics seeds; separate full-window observation MC diagnostic is available.',
            'Nested thinning uses correlated samples; independent pairs remain different.',
            '4-Gyr preflight is sparse, fixed-radius and conditional; zero sampled contribution does not prove physical zero.',
            'Exact-template profile calibration excludes dynamics/population/background misspecification and does not guarantee uniform nuisance coverage.',
            'Excluded duplicate-seed startup attempt is not used; follow-up explicitly uses seed271828.'])
    (out/'summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    import os
    os.environ.setdefault('MPLCONFIGDIR',str(root/'.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    axes[0].errorbar(ages.T_myr/1000,ages['mean'],yerr=ages.SE,fmt='o-')
    axes[0].set(xlabel='Lookback window (Gyr)',ylabel='Expected selected stars',title='Sparse age preflight, prescribed fixed radius')
    x=np.arange(len(enlarged));means=np.array([r['coverage'] for r in enlarged])
    lower=np.array([r['wilson_95'][0] for r in enlarged]);upper=np.array([r['wilson_95'][1] for r in enlarged])
    axes[1].errorbar(x,means,yerr=[means-lower,upper-means],fmt='o',capsize=4)
    axes[1].axhline(.95,color='black',ls='--');axes[1].set_xticks(x,['Rate','Time','Time + rate','Full point'])
    axes[1].set(ylabel='Fresh held-out inclusion fraction',title='1000 calibration + 500 new tests',ylim=(.85,1.))
    fig.suptitle('Conditional numerical experiments, not a validated physical mass-loss measurement',fontsize=11)
    fig.savefig(out/'followup_diagnostics.png',dpi=160);plt.close(fig)
    lines=['# 改善方法的實際測試','',
        '這次完成三項測試：三個獨立種子各 8,000 粒子的空間收斂；延伸至 4 Gyr 的稀疏年代預演；分開統計量並擴大校準的 profile 回收。尚未宣稱絕對物理 MLR 已驗證。','',
        '## 8,000 粒子','',
        '本次固定 MLR = 9 M☉/Myr、人口 tilt = 0，T = 250、500、750、1000 Myr。原細像素約 0.111°，四個 T 中兩個通過嚴格主要像素 MC 目標；約 0.444° 的像素四個都通過。不能把單一流失率的結果直接推廣到其他 MLR 或星團。','',
        '| 像素 deg | T Myr | 最差模板 SE/Poisson | 通過 |', '|---|---|---|---|']
    for r in noise.loc[noise.factor.isin([1,4])].to_dict('records'):
        lines.append(f"| {r['cell_deg']:.3f} | {r['T_myr']:g} | {r['maximum_se_over_poisson']:.3f} | {r['major_pixel_target_passed']} |")
    lines+=['',
        '4,000→8,000 的巢狀 TV 約 0.044–0.057，但這些層級相關。真正獨立 8,000 粒子種子配對的 TV 約 0.110–0.148，仍不全低於 0.1。兩者不能混稱為獨立收斂證明。',
        '在固定 rate 與人口的時間 profile 中，巢狀 4k→8k 的最大 q 改變，在細網格約 2.50、0.444° 網格約 0.13；更粗分箱也會改變資料資訊量，不能只為了通過而合併。',
        '', '[完整密度診斷](../Pal5_dense_8000_analysis_v1/REPORT_ZH.md)','',
        '## 年代範圍','',
        '每模型僅 400 粒子、三種子；僅用於年代貢獻預演。相同指定質量歷史下，1 Gyr 預期約 79.86 顆，1.5 Gyr 約 85.15 顆，2–4 Gyr 的取樣模型沒有再增加。1 Gyr 外人口約占此 4 Gyr 模型訊號 6.2%。',
        '建議下一階段先以高取樣確認 1.5–2 Gyr。零取樣貢獻不證明真實貢獻精確為零；4 Gyr 過去質量約為現在 3.69 倍，固定尺度假設也需額外驗證。','',
        '## Profile 校準','',
        '原先完整候選點、MLR profile、T profile 與 T–MLR 聯合區域已分開計算。先作三情境各 200 次校準、300 次測試，再將舊 T500 的 500 個樣本改作校準／開發資料，加入獨立種子 271828 的 500 個校準及 500 個全新測試樣本。最終門檻使用 1,000 校準分數的第 951 小值，沒有依全新測試結果調整。','',
        '| 統計目標 | 全新 500 次涵蓋率 | Wilson 95% 區間 |','|---|---|---|']
    for row in enlarged:
        lines.append(f"| {row['target']} | {row['coverage']:.1%} | {row['wilson_95'][0]:.1%}–{row['wilson_95'][1]:.1%} |")
    lines+=['',
        '這些結果與名義 95% 相容，但只是 rate9、T500、指定人口與背景下的精確模板測試；沒有驗證全部干擾參數、非網格真值、動力學失配或完整選星重建。Wilson 區間描述固定所得門檻的測試比例，不包含全部校準批次不確定性。','',
        '## 下一步','',
        '以 8,000 粒子／種子、三種子及約 0.444° 的總量分析作下一輪候選設定，對所有候選 MLR 一致測試；另保留細網格作形態診斷。先高取樣確認至 1.5–2 Gyr，再加密 MLR 網格。不能只將最佳 MLR 換成高取樣模板、其他候選留在低取樣，然後直接比較。',
        '細尾巴寬度、人口質量換算、背景與半徑演化仍需獨立驗證。60 項測試通過不等於這些物理假設已通過。','',
        '![Age and calibration](followup_diagnostics.png)','',
        '[機器可讀摘要](summary.json) · [測試紀錄](unittest.log)']
    (out/'REPORT_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    (out/'completion.json').write_text(json.dumps(dict(completed=True,physical_rate_measurement_validated=False,
        report_sha256=sha256_file(out/'REPORT_ZH.md'),summary_sha256=sha256_file(out/'summary.json')),indent=2))
    print('Improvement experiments summarized:',out)
    return result


if __name__=='__main__':run()
