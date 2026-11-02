"""Four MLR comparisons for a selected sample, optionally including NGC362.

Keep the primary integrated estimates and mark NGC362's score-only estimate.
Also export a uniform score-profile version without changing the main results.
"""
import argparse
import json
import math
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from plot_mlr_parameters import add_labels, q_regions, read, sha

IDS=('NGC288','NGC5466','NGC6341','NGC7099','Pal5','NGC5904','NGC7078','NGC3201','NGC362')
AXES=(('mass_1e5_msun','星團質量',r'質量 $M$ ($10^5\,M_\odot$)','mass'),
      ('specific_rate_gyr_inv','MLR／M',r'MLR / $M$ ($\mathrm{Gyr}^{-1}$)','specific_rate'),
      ('R_gc_kpc','銀心距離',r'銀心距離 $R_{\rm GC}$ (kpc)','galactocentric_distance'),
      ('eccentricity','軌道偏心率',r'軌道偏心率 $e$','eccentricity'))
BLUE='#286FA0'
ORANGE='#C27023'
PURPLE='#86559C'


def snapshot_source(path,digest,root,manifest):
    """Use the copied snapshot first, so live outputs cannot change a figure."""
    record=next((r for r in manifest if r['source']==path and r['sha256']==digest),None)
    candidate=root/record['destination'] if record is not None else Path(path)
    if not candidate.is_file() or sha(candidate)!=digest:raise ValueError('Source snapshot identity failed: '+path)
    return candidate


def derive(row,master,profile,method):
    row=dict(row);cid=row['cluster'];cat=master.loc[cid]
    x=profile.rate_msun_myr.to_numpy();q=profile.q.to_numpy();regions=q_regions(x,q)
    best=float(row['rate']);chosen=next(r for r in regions if r['lower']-1e-9<=best<=r['upper']+1e-9)
    if abs(q[np.flatnonzero(x==best)[0]])>1e-7:raise ValueError('Best rate is not a profile minimum')
    mass=float(cat['mass']);ecc=float((cat.R_apo-cat.R_peri)/(cat.R_apo+cat.R_peri))
    if not math.isclose(ecc,float(cat.ecc),abs_tol=5.1e-6):raise ValueError('Eccentricity definition differs')
    row.update(label='NGC '+cid[3:] if cid.startswith('NGC') else 'Pal 5',method=method,delta=0.,
        rate_msun_myr=best,mass_msun=mass,mass_1e5_msun=mass/1e5,R_gc_kpc=float(cat.R_gc),
        eccentricity=float(cat.ecc),specific_rate_gyr_inv=1000*best/mass,
        q1_regions=regions,q1_lower=chosen['lower'],q1_upper=chosen['upper'],
        diagnostic_minus=None if chosen['lower_grid_truncated'] else best-chosen['lower'],
        diagnostic_plus=None if chosen['upper_grid_truncated'] else chosen['upper']-best,
        rate_best_on_grid_edge=best in (x[0],x[-1]))
    return row


def ols(rows,xkey):
    x=np.array([r[xkey] for r in rows]);y=np.array([r['rate_msun_myr'] for r in rows])
    a=float(np.sum((x-x.mean())*(y-y.mean()))/np.sum((x-x.mean())**2));b=float(y.mean()-a*x.mean())
    check=np.linalg.lstsq(np.column_stack([x,np.ones(len(x))]),y,rcond=None)[0]
    if not np.allclose([a,b],check,atol=1e-10,rtol=1e-10):raise ValueError('Independent OLS formula mismatch')
    residual=float(np.sum((y-a*x-b)**2));total=float(np.sum((y-y.mean())**2))
    r2=1-residual/total
    if not np.isclose(r2,np.corrcoef(x,y)[0,1]**2,atol=1e-10):raise ValueError('R squared cross-check failed')
    return dict(x_variable=xkey,a=a,b=b,r_squared=r2,n=len(rows),
        fit_clusters=[r['cluster'] for r in rows],weighting='unweighted ordinary least squares',
        includes_grid_boundary_estimates=any(r['rate_best_on_grid_edge'] for r in rows),
        mathematically_coupled_axes=xkey=='specific_rate_gyr_inv',confidence_intervals_calibrated=False)


def panel(ax,rows,xkey,title,xlabel,fit):
    diag=[]
    for r in rows:
        for region in r['q1_regions']:
            diag.append((r,region))
    if xkey=='specific_rate_gyr_inv':
        xmax=max(1000*reg['upper']/r['mass_msun'] for r,reg in diag);ax.set_xlim(0,xmax*1.12)
    elif xkey=='mass_1e5_msun':ax.set_xlim(0,max(r[xkey] for r in rows)*1.13)
    else:
        xmin=min(r[xkey] for r in rows);xmax=max(r[xkey] for r in rows);margin=(xmax-xmin)*.12
        ax.set_xlim(max(0,xmin-margin),min(1,xmax+margin) if xkey=='eccentricity' else xmax+margin)
    ymax=max(reg['upper'] for _,reg in diag);ax.set_ylim(0,ymax*1.17)
    ax.set(title=title,xlabel=xlabel,ylabel=r'MLR ($M_\odot\,\mathrm{Myr}^{-1}$)')
    ax.grid(True,alpha=.18,linewidth=.7);ax.set_axisbelow(True)
    for side in ['top','right']:ax.spines[side].set_visible(False)
    for r in rows:
        special=r['cluster']=='NGC362';highlight=r['cluster']=='NGC6139'
        color=ORANGE if special else PURPLE if highlight else BLUE
        for region in r['q1_regions']:
            lo,hi=region['lower'],region['upper']
            if xkey=='specific_rate_gyr_inv':
                # Both axes contain the SAME MLR. The allowed pairs form a line,
                # not independent horizontal and vertical error bars.
                xx=[1000*lo/r['mass_msun'],1000*hi/r['mass_msun']]
            else:xx=[r[xkey],r[xkey]]
            ax.plot(xx,[lo,hi],color=color,lw=1.2,zorder=2)
            for xend,yend,truncated in [(xx[0],lo,region['lower_grid_truncated']),(xx[1],hi,region['upper_grid_truncated'])]:
                if not truncated:ax.plot(xend,yend,marker='_',color=color,ms=6,mew=1.2,zorder=2)
        ax.plot(r[xkey],r['rate_msun_myr'],marker='D' if special else 's' if highlight else 'o',ms=7.3 if special or highlight else 6.4,lw=0,
            markerfacecolor='white' if r['rate_best_on_grid_edge'] else color,markeredgecolor=color,markeredgewidth=1.3,zorder=4)
    extra=[(dx,dy) for dx in (60,90,120) for dy in (0,20,-20,40,-40,60)]
    add_labels(ax,rows,xkey,extra_offsets=extra,legacy_offsets=False)
    fit_rows=[r for r in rows if r['cluster'] in fit['fit_clusters']]
    xx=np.linspace(min(r[xkey] for r in fit_rows),max(r[xkey] for r in fit_rows),200)
    ax.plot(xx,fit['a']*xx+fit['b'],ls='--',color='#59616C',lw=1.25,zorder=1)
    sign='+' if fit['b']>=0 else '-'
    text=rf"$y={fit['a']:.3f}x{sign}{abs(fit['b']):.3f}$"+'\n'+rf"$R^2={fit['r_squared']:.3f}$"+f"，擬合 N={fit['n']}"
    ax.text(.97,.97,text,transform=ax.transAxes,ha='right',va='top',fontsize=10,
        bbox=dict(facecolor='white',edgecolor='#D7DDE3',boxstyle='round,pad=.4',alpha=.96))


def run(snapshot,ngc362,out,clusters=IDS):
    snapshot,ngc362,out=map(lambda p:Path(p).resolve(),(snapshot,ngc362,out))
    ids=tuple(clusters)
    if len(ids)<3 or len(set(ids))!=len(ids):raise ValueError('Choose at least three distinct clusters')
    if out.exists():raise FileExistsError('Choose a fresh plot version')
    data=read(snapshot/'data.json');manifest=read(snapshot/'source_manifest.json')['files']
    master_path=snapshot/'catalogue/master_clusters.csv';master=pd.read_csv(master_path).set_index('id')
    catalog=next(r for r in manifest if r['destination']=='catalogue/master_clusters.csv')
    if sha(master_path)!=catalog['sha256']:raise ValueError('Frozen physical catalogue changed')
    primary=[];raw=[];sources=[]
    for cid in ids:
        if cid=='NGC362':
            r=next(r for r in read(ngc362/'data.json')['rows'] if r['delta']==0.)
            original=Path(r['source_json']);source=ngc362/'source_scores/offset 0'/original.name
            if sha(source)!=r['source_sha256']:raise ValueError('NGC362 source score changed')
            base=dict(cluster=cid,rate=r['rate']);digest=r['source_sha256']
        else:
            r=next(r for r in data['rows'] if r['cluster']==cid and r['delta']==0.)
            profile=pd.read_csv(snapshot/'completed'/cid/'offset_0/rate_profile.csv')
            primary.append(derive(r,master,profile,'background_and_MF_integrated'))
            original,digest=next((p,h) for p,h in r['source_hashes'].items() if p.endswith('joint_grid.json'))
            source=snapshot_source(original,digest,snapshot,manifest);base=dict(cluster=cid)
        score=read(source);frame=pd.DataFrame(score['candidates'])
        if len(frame)!=1600 or frame.duplicated(['rate_msun_myr','duration_myr','scenario']).any():raise ValueError('Incomplete score grid')
        if any(json.loads(s)['escape_slope_offset']!=0 for s in frame.scenario):raise ValueError('Nonzero offset in plot')
        if not np.allclose(frame.q,2*(frame.negative_log_likelihood-frame.negative_log_likelihood.min()),atol=1e-9):raise ValueError('Raw q mismatch')
        best=frame.loc[frame.negative_log_likelihood.idxmin()];base['rate']=float(best.rate_msun_myr)
        profile=frame.groupby('rate_msun_myr',as_index=False).q.min()
        raw_row=derive(base,master,profile,'profile_over_T_and_five_population_nodes_background_optimized')
        raw.append(raw_row)
        if cid=='NGC362':primary.append(raw_row)
        sources.append(dict(cluster=cid,score_snapshot=str(source),score_sha256=digest))
    if [r['cluster'] for r in primary]!=list(ids) or len(raw)!=len(ids):raise ValueError('Requested cluster membership differs')
    out.mkdir(parents=True);os.environ['MPLCONFIGDIR']=str(out/'mplcache')
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.lines import Line2D
    font_manager.fontManager.addfont('C:/Windows/Fonts/msjh.ttc')
    plt.rcParams.update({'font.family':'Microsoft JhengHei','axes.unicode_minus':False,'font.size':10,
        'axes.titlesize':12,'axes.labelsize':11,'pdf.fonttype':42,'ps.fonttype':42})
    reports={}
    for name,rows,uniform in [('current_estimates',primary,False),('uniform_raw_scores',raw,True)]:
        dest=out/name;dest.mkdir()
        fit_rows=rows if uniform else [r for r in rows if r['cluster']!='NGC362']
        fits={key:ols(fit_rows,key) for key,_,_,_ in AXES}
        pd.DataFrame(fits.values()).to_csv(dest/'linear_fits.csv',index=False,encoding='utf-8-sig')
        pd.DataFrame([{k:v for k,v in r.items() if k not in ('q1_regions','source_hashes','q_regions')} for r in rows]).to_csv(dest/'plot_data.csv',index=False,encoding='utf-8-sig')
        note=(r'診斷線段：$q_{\rm prof}\leq1$，網格間線性插值；不是標準差或信賴區間。'
            +'\nMLR／M 圖的斜線段表示同一 MLR 在兩軸同步變動；未傳播質量誤差。'
            +'\nMLR／M 與 MLR 數學上耦合，其相關性不能當作獨立物理證據。')
        main_count=len(rows)-int('NGC362' in ids)-int('NGC6139' in ids)
        note+=(f'\n虛線：{len(rows)} 顆同定義原始 scores 的未加權 y=ax+b。' if uniform else f'\n虛線：只擬合 {len(fit_rows)} 顆積分比較；NGC 362 的原始 scores 點未混入擬合。' if 'NGC362' in ids else f'\n虛線：{len(fit_rows)} 顆積分比較的未加權 y=ax+b。')
        legend=[Line2D([],[],marker='o',color=BLUE,lw=0,label=f'其他 {main_count} 顆（原始 scores）' if uniform else f'{main_count} 顆（積分比較）')]
        if 'NGC6139' in ids:legend.append(Line2D([],[],marker='s',color=PURPLE,lw=0,label='NGC 6139（原始 scores）' if uniform else 'NGC 6139（積分比較）'))
        if 'NGC362' in ids:legend.append(Line2D([],[],marker='D',color=ORANGE,lw=0,label='NGC 362（原始 scores）'))
        legend.append(Line2D([],[],marker='o',mfc='white',mec='#59616C',lw=0,label='網格邊界最佳點；無端帽＝截斷'))
        for key,title,xlabel,filename in AXES:
            fig,ax=plt.subplots(figsize=(9,7.3));fig.subplots_adjust(left=.11,right=.98,bottom=.26,top=.82)
            panel(ax,rows,key,'MLR 與'+title,xlabel,fits[key])
            fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.54,.975),ncol=2,frameon=False,fontsize=9)
            fig.text(.98,.875,f'offset＝0；繪圖 N＝{len(rows)}',ha='right',fontsize=10,color='#53616F')
            fig.text(.11,.025,note,fontsize=9,color='#53616F',linespacing=1.6,va='bottom')
            for ext in ('png','pdf','svg'):fig.savefig(dest/('mlr_vs_'+filename+'.'+ext),dpi=220)
            plt.close(fig)
        fig,axes=plt.subplots(2,2,figsize=(15.8,12.2));fig.subplots_adjust(left=.075,right=.98,bottom=.17,top=.89,wspace=.23,hspace=.36)
        for ax,(key,title,xlabel,_) in zip(axes.flat,AXES):panel(ax,rows,key,'MLR 與'+title,xlabel,fits[key])
        fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.53,.985),ncol=2 if len(legend)>3 else 3,frameon=False,fontsize=10)
        title=f'{len(rows)} 顆均採原始 scores 的剖面比較' if uniform else f'保留目前估計值：{len(fit_rows)} 顆積分比較'+('＋NGC 362 原始 scores' if 'NGC362' in ids else '')
        fig.text(.53,.94,title,ha='center',fontsize=12,color='#334155')
        fig.text(.075,.025,note,fontsize=10,color='#53616F',linespacing=1.6,va='bottom')
        for ext in ('png','pdf','svg'):fig.savefig(dest/('four_parameter_comparison.'+ext),dpi=190)
        plt.close(fig)
        reports[name]=dict(records=rows,fits=fits,plot_clusters=len(rows),fit_clusters=len(fit_rows))
    report=dict(created_taipei=datetime.now().astimezone().isoformat(),requested_clusters=list(ids),delta=0.,
        source_catalogue=str(master_path),source_catalogue_sha256=sha(master_path),source_scores=sources,variants=reports,
        uncertainty_definition='q<=1 profile intervals, not standard deviations or confidence intervals',
        ratio_units='MLR/M in Gyr^-1, from 1000*MLR/(adopted present-day catalogue mass)',
        ratio_uncertainty='Both axes use the same MLR; diagonal profile segments hold catalogue M fixed. Mass uncertainty is not propagated.',
        ratio_mass_scope='Escaping represented stellar population and adopted catalogue total mass do not have identical mass scopes.',
        target_posterior_used=False,confidence_intervals_calibrated=False,physical_rate_measurement_validated=False)
    (out/'data_and_provenance.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    text=[f'# 指定 {len(ids)} 顆星團的 MLR 比較','',
        '星團：'+'、'.join(ids)+'。全為 offset=0。','',
        '縱軸皆為 MLR（M☉/Myr）；四橫軸依序為 M/(10^5 M☉)、MLR/M（Gyr⁻¹）、目前三維銀心距離（kpc）、軌道偏心率。','',
        '## 兩個版本','',
        f'- current_estimates/：沿用既有積分比較的 MLR；若納入 NGC362，其原始 scores 用橘色菱形另標。繪圖 N={len(ids)}，虛線只對同方法的{len(primary)-int("NGC362" in ids)}顆作未加權 y=ax+b。',
        f'- uniform_raw_scores/：{len(ids)}顆全部重取原始 scores 的最佳值，背景先最佳化，再對 T 與五個人口節點取最大概似，虛線以{len(ids)}顆擬合。此版本改變估計定義，不取代原有積分結果；例如 Pal5 的最佳 MLR 為9，主版本為5。','',
        '每版均包含四張獨立圖及四圖合併，提供 PNG、PDF、SVG；plot_data.csv 與 linear_fits.csv 保留完整精度。','',
        '## 解讀','',
        'q≤1 棒為診斷範圍，不是標準差或經校準的信賴區間。空心點代表最佳值在計算網格邊界；無端帽的方向被網格截斷，不表示物理上下限。',
        'MLR/M 圖的斜線段沿 (MLR/M,MLR) 同步變化，避免把同源的兩個量畫成獨立的橫向與縱向誤差。質量採固定目錄值，未傳播質量誤差。',
        'MLR 與 MLR/M 具有數學耦合，因此此圖的 R² 不能作獨立物理關係的證據。所代表逃逸人口質量範圍與目錄總質量也不完全一致，比值為條件式比較指標。',
        '擬合是描述性普通最小平方，包含網格邊界最佳值，不以 q 棒加權，不提供顯著性、因果或正式信賴區間。',
        '銀心距離、質量與軌道同時參與模型設定，不能將相關性視為獨立模型驗證。僅繪製縱向 q 診斷不代表橫軸參數沒有誤差。',
        '保留各來源 SHA256；未重新跑動力學，未重新下載或更換目錄值。']
    (out/'README.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    print(json.dumps(dict(out=str(out),requested_clusters=ids,panels=8,variants={k:dict(plot=v['plot_clusters'],fit=v['fit_clusters']) for k,v in reports.items()}),ensure_ascii=True))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--snapshot',type=Path,default=Path('outputs/consolidated_results_20261102_v2'))
    p.add_argument('--ngc362',type=Path,default=Path('outputs/NGC362_score_summary_20261102_v1'))
    p.add_argument('--clusters',nargs='+',default=list(IDS))
    p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.snapshot,a.ngc362,a.out,a.clusters)
