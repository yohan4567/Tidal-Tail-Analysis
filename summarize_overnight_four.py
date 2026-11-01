"""Read-only summary of the four November 1–2 completed model comparisons."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOTS={
    'M5 (NGC 5904)':HERE/'outputs/s7r2/NGC5904/local_v3/results',
    'NGC 3201':HERE/'outputs/b3201/results',
    'M92 (NGC 6341)':HERE/'outputs/q7r/NGC6341/results',
    'M30 (NGC 7099)':HERE/'outputs/q7r/NGC7099/results',
}


def run():
    out=HERE/'outputs/overnight_four_20261102';out.mkdir(exist_ok=True)
    records=[]
    for cluster,root in ROOTS.items():
        if json.loads((root/'run_status.json').read_text())['stage']!='completed_model_comparisons_validation_pending':
            raise ValueError('Incomplete cluster: '+cluster)
        for offset,delta in [('offset_0',0.),('offset_minus0p5',-.5),('offset_minus1',-1.)]:
            path=root/offset/'comparison/comparison.json';data=json.loads(path.read_text());node=data['best_sampled_node']
            records.append(dict(cluster=cluster,delta=delta,**node,lookback_gyr=node['duration_myr']/1000,
                rate_on_grid_boundary=node['rate_msun_myr'] in (1.,31.),time_on_grid_boundary=node['duration_myr'] in (100.,2000.),
                result_folder=str(root),source_json=str(path),confidence_intervals_calibrated=False,physical_rate_measurement_validated=False))
    pd.DataFrame(records).to_csv(out/'summary.csv',index=False,encoding='utf-8-sig')
    (out/'summary.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for ax,(cluster,root) in zip(axes.flat,ROOTS.items()):
        data=pd.read_csv(root/'offset_0/comparison/likelihood_grid.csv')
        grid=data.pivot(index='rate_msun_myr',columns='duration_myr',values='q').sort_index()
        plot=ax.pcolormesh(grid.columns.to_numpy()/1000,grid.index.to_numpy(),grid.to_numpy(),shading='nearest',
                           cmap='viridis_r',vmin=0,vmax=20,rasterized=True)
        best=next(r for r in records if r['cluster']==cluster and r['delta']==0.)
        ax.scatter(best['lookback_gyr'],best['rate_msun_myr'],marker='*',s=130,facecolor='white',edgecolor='black',linewidth=.8)
        ax.set(title=cluster,xlabel='Lookback window T (Gyr)',ylabel=r'Escape rate ($M_\odot$ / Myr)')
        ax.text(.97,.96,f"Best grid point: R = {best['rate_msun_myr']:g}, T = {best['lookback_gyr']:g} Gyr",
                transform=ax.transAxes,ha='right',va='top',fontsize=9,bbox=dict(facecolor='white',alpha=.9,edgecolor='none',pad=3))
    cb=fig.colorbar(plot,ax=axes,extend='max',shrink=.9);cb.set_label(r'$q=2(\ell_{\max}-\ell)$; lower is better')
    fig.suptitle(r'Completed conditional model comparisons: $\delta=0$'+'\nGrid best fits; no calibrated confidence intervals',fontsize=12)
    fig.savefig(out/'offset0_comparison.png',dpi=190);fig.savefig(out/'offset0_comparison.pdf');plt.close(fig)
    lines=['# 昨晚完成的四顆星團','',
           '數值為對背景與質量函數量測不確定性積分後，在既有網格中選出的最佳模型。未計算 R/T 後驗，亦無經校準的信賴區間。',
           'R 單位：M☉/Myr；T 單位：Gyr。每個情境有 320 個 R–T 網格點。','',
           '| 星團 | δ | R | T | 網格邊界 |','|---|---:|---:|---:|---|']
    for r in records:
        boundary=[]
        if r['rate_on_grid_boundary']:boundary.append('R 下界' if r['rate_msun_myr']==1 else 'R 上界')
        if r['time_on_grid_boundary']:boundary.append('T 下界' if r['duration_myr']==100 else 'T 上界')
        lines.append(f"| {r['cluster']} | {r['delta']:g} | {r['rate_msun_myr']:g} | {r['lookback_gyr']:g} | {'、'.join(boundary) or '未落在邊界'} |")
    lines+=['','M5 的 δ=0、−0.5 最佳 R 落在網格下界 1，現有網格不足以定位較低流失率；不得當作精確測量或以零誤差呈現。',
            '圖中的 q 僅為模型相對分數，色階在 20 以上飽和，不代表 95% 信賴區域。T 是逃逸歷史的回溯視窗，不是星團年齡。','',
            '![δ=0 模型比較](offset0_comparison.png)','','## 原始結果資料夾','']
    for cluster,root in ROOTS.items():lines.append(f'- [{cluster}](<{root.as_posix()}>)')
    (out/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(out)


if __name__=='__main__':run()
