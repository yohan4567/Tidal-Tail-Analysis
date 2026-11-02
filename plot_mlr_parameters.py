"""Plot delta=0 integrated-likelihood grid estimates with q<=1 diagnostic spans.

No posterior or confidence interval is computed. Piecewise-linear profile
interpolation is a display approximation, not a finer dynamical model grid.
"""
import argparse,hashlib,json,math,os
from pathlib import Path
from datetime import datetime,timezone,timedelta
import zipfile,xml.etree.ElementTree as ET
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
STUDY=HERE.parent/'TidaleffectsOnGCs-main/TidaleffectsOnGCs-main/tidal_study_v2'
LABELS={'NGC288':'NGC 288','Pal5':'Pal 5','NGC6171':'M107','NGC5466':'NGC 5466',
    'NGC6144':'NGC 6144','NGC2298':'NGC 2298','NGC5904':'M5','NGC3201':'NGC 3201',
    'NGC6341':'M92','NGC7099':'M30','NGC7078':'M15','NGC6121':'M4','NGC6397':'NGC 6397',
    'NGC7089':'M2','NGC6139':'NGC 6139','NGC6266':'M62','NGC362':'NGC 362'}
COMPARISONS={
    **{c:HERE/'outputs/priority_six_20261101_v1'/c/'results/offset_0/comparison' for c in ['NGC5466','NGC6144','NGC2298']},
    'NGC5904':HERE/'outputs/s7r2/NGC5904/local_v3/results/offset_0/comparison',
    'NGC3201':HERE/'outputs/b3201/results/offset_0/comparison',
    **{c:HERE/'outputs/q7r'/c/'results/offset_0/comparison' for c in ['NGC6341','NGC7099','NGC7078','NGC6139','NGC6266']},
    **{f'NGC{c}':HERE.parent/'results'/c/'0/comparison' for c in ['6121','6397','7089']},
    'NGC362':HERE.parent/'results/ngc362/offset 0/comparison',
}

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def canonical(value):return value.replace(' ','').replace('_','')

def document_tables(path):
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    with zipfile.ZipFile(path) as z:root=ET.fromstring(z.read('word/document.xml'))
    tables=[]
    for table in root.findall('.//w:body/w:tbl',ns):
        tables.append([[''.join(n.text or '' for n in c.findall('.//w:t',ns)) for c in r.findall('w:tc',ns)] for r in table.findall('w:tr',ns)])
    shapes={};physical={}
    for row in tables[0][1:]:
        for offset in (0,4):
            if row[offset].strip():shapes[canonical(row[offset])]=dict(S1=float(row[offset+1]),S2=float(row[offset+2]),shape=row[offset+3])
    for row in tables[1][1:]:
        for offset in (0,5):
            if row[offset].strip():physical[canonical(row[offset])]=dict(mass=float(row[offset+1].replace(',','')),
                ecc=float(row[offset+2]),R_gc=float(row[offset+4]))
    if len(shapes)!=29 or len(physical)!=29:raise ValueError('Expected 29 unique cluster rows in both supplied tables')
    return shapes,physical

def q_regions(rates,q,threshold=1.):
    """All connected regions of a piecewise-linear diagnostic profile."""
    x=np.asarray(rates,float);y=np.asarray(q,float)
    if len(x)<2 or not np.isfinite([x,y]).all() or (np.diff(x)<=0).any() or (y< -1e-8).any():
        raise ValueError('Invalid sorted diagnostic profile')
    intervals=[]
    for i in range(len(x)-1):
        a,b,qa,qb=x[i],x[i+1],y[i],y[i+1]
        if qa<=threshold and qb<=threshold:segment=[a,b]
        elif qa>threshold and qb>threshold:continue
        else:
            cross=a+(threshold-qa)*(b-a)/(qb-qa)
            segment=[a,cross] if qa<=threshold else [cross,b]
        if intervals and math.isclose(intervals[-1][1],segment[0],abs_tol=1e-10):intervals[-1][1]=segment[1]
        else:intervals.append(segment)
    return [dict(lower=float(a),upper=float(b),lower_grid_truncated=bool(a==x[0] and y[0]<=threshold),
                 upper_grid_truncated=bool(b==x[-1] and y[-1]<=threshold)) for a,b in intervals]

def load_estimates():
    estimates={};excluded={}
    old=read(HERE/'outputs/offset_results_20261101/summary.json')
    for r in old['results']:
        if r['delta']!=0 or r['cluster'] not in ('NGC288','Pal5','NGC6171'):continue
        for p,digest in r['sources'].items():
            if sha(p)!=digest:raise ValueError('Legacy integrated source changed: '+p)
        estimates[r['cluster']]=dict(best=r['best'],profile=r['rate_profile'],sources=r['sources'])
    for cid,folder in COMPARISONS.items():
        report=folder/'comparison.json'
        if not report.exists():excluded[cid]='尚未取得 offset=0 的干擾參數積分比較結果';continue
        meta=read(report)
        if meta['plan']['method']!='nuisance_integrated_likelihood_no_target_posterior_v1' or meta['target_parameter_posterior_computed']:
            raise ValueError('Incompatible result definition')
        grid=pd.read_csv(folder/'likelihood_grid.csv');profile=pd.read_csv(folder/'rate_profile.csv').sort_values('rate_msun_myr')
        q=2*(grid.log_likelihood.max()-grid.log_likelihood)
        if len(grid)!=320 or grid.duplicated(['rate_msun_myr','duration_myr']).any() or not np.allclose(grid.q,q,atol=1e-9):raise ValueError('Invalid comparison grid')
        expected=grid.groupby('rate_msun_myr',as_index=False).q.min().sort_values('rate_msun_myr')
        pd.testing.assert_frame_equal(profile.reset_index(drop=True),expected.reset_index(drop=True),check_dtype=False,atol=1e-9,rtol=1e-10)
        source=next((p for p in [folder.parent/'grid/scores/joint_grid.json',folder.parent/'scores/joint_grid.json'] if p.exists()),None)
        if source is None or sha(source)!=meta['plan']['score_sha256']:raise ValueError('Comparison source score identity mismatch')
        if any(json.loads(r['scenario'])['escape_slope_offset']!=0 for r in read(source)['candidates']):raise ValueError('Result is not offset=0')
        if meta['plan']['beta_weight']!={'distribution':'uniform','lower':0.,'upper':5.}:raise ValueError('Background integration range differs')
        estimates[cid]=dict(best=meta['best_sampled_node'],profile=profile.to_dict('records'),
                           sources={str(p):sha(p) for p in [report,folder/'likelihood_grid.csv',folder/'rate_profile.csv',source]})
    return estimates,excluded

def add_labels(ax,rows,xkey):
    """Choose readable deterministic offsets; keep every plotted cluster label."""
    from matplotlib.transforms import Bbox
    ax.figure.canvas.draw();renderer=ax.figure.canvas.get_renderer();occupied=[]
    points=[ax.transData.transform((r[xkey],r['rate_msun_myr'])) for r in rows]
    offsets=[(7,9),(7,-15),(-7,9),(-7,-15),(0,22),(0,-27),(23,3),(-23,3),
             (22,20),(-22,20),(22,-25),(-22,-25),(0,38),(0,-42),(38,10),(-38,10)]
    for r in sorted(rows,key=lambda r:(r['rate_msun_myr'],r[xkey])):
        point=(r[xkey],r['rate_msun_myr']);choice=None
        candidates=offsets
        if xkey=='R_gc_kpc' and r['cluster']=='NGC6121':candidates=[(8,15)]
        if xkey=='R_gc_kpc' and r['cluster']=='NGC5904':candidates=[(-8,-18)]
        for dx,dy in candidates:
            ha='left' if dx>0 else ('right' if dx<0 else 'center')
            text=ax.annotate(r['label'],point,xytext=(dx,dy),textcoords='offset points',ha=ha,va='center',fontsize=9)
            box=text.get_window_extent(renderer).expanded(1.18,1.25);text.remove()
            area=sum(max(0,min(box.x1,b.x1)-max(box.x0,b.x0))*max(0,min(box.y1,b.y1)-max(box.y0,b.y0)) for b in occupied)
            npoints=sum(box.contains(*p) for p in points)
            outside=max(0,ax.bbox.x0-box.x0)+max(0,box.x1-ax.bbox.x1)+max(0,ax.bbox.y0-box.y0)+max(0,box.y1-ax.bbox.y1)
            cost=area*100+npoints*5000+outside*10000+abs(dx)+abs(dy)
            if choice is None or cost<choice[0]:choice=(cost,dx,dy,ha,box)
        _,dx,dy,ha,box=choice;occupied.append(box)
        ax.annotate(r['label'],point,xytext=(dx,dy),textcoords='offset points',ha=ha,va='center',fontsize=9,
            bbox=dict(facecolor='white',edgecolor='none',alpha=.88,pad=.5),
            arrowprops=dict(arrowstyle='-',color='#7E8790',lw=.6,shrinkA=2,shrinkB=5))

def linear_fit(rows,key):
    x=np.array([r[key] for r in rows]);y=np.array([r['rate_msun_myr'] for r in rows])
    if np.unique(x).size<2:raise ValueError('Linear fit requires varying x values')
    a,b=np.linalg.lstsq(np.column_stack([x,np.ones(len(x))]),y,rcond=None)[0]
    residual=float(np.sum((y-(a*x+b))**2));total=float(np.sum((y-y.mean())**2))
    # Independently check the closed-form OLS slope and intercept.
    expected=np.sum((x-x.mean())*(y-y.mean()))/np.sum((x-x.mean())**2)
    if not np.isclose(a,expected,rtol=1e-10,atol=1e-12) or not np.isclose(b,y.mean()-a*x.mean(),atol=1e-10):raise ValueError('OLS cross-check failed')
    return dict(x_variable=key,a=float(a),b=float(b),r_squared=1-residual/total if total>0 else None,
        residual_sum_squares=residual,n=len(rows),weighting='unweighted ordinary least squares',
        grid_boundary_best_points=sum(r['rate_best_on_grid_edge'] for r in rows),
        interpretation='Descriptive fit to sampled conditional estimates, including grid-edge best points; q<=1 spans are not standard deviations and are not used as weights.')

def panel(ax,rows,xkey,title,xlabel,fit=None):
    colors={'雙側':'#2878AD','S 形':'#B87520','單側':'#8E5196'}
    markers={'雙側':'o','S 形':'s','單側':'^'}
    if xkey=='mass_msun':ax.set_xscale('log');ax.set_xlim(9000,850000)
    elif xkey=='mass_1e5_msun':ax.set_xlim(0,7.)
    elif xkey=='R_gc_kpc':ax.set_xlim(0,19.5)
    else:ax.set_xlim(.12,1.)
    ymax=max(r['q1_upper'] for r in rows);ax.set_ylim(.35,max(9.,ymax+1.))
    ax.set(title=title,xlabel=xlabel,ylabel=r'MLR ($M_\odot\,\mathrm{Myr}^{-1}$)')
    ax.grid(True,alpha=.18,lw=.7);ax.set_axisbelow(True)
    for side in ['top','right']:ax.spines[side].set_visible(False)
    for r in rows:
        color=colors[r['tail_shape']];x=r[xkey]
        # Uncapped endpoints indicate an unresolved grid-truncated q region.
        for region in r['q1_regions']:
            lo,hi=region['lower'],region['upper'];ax.plot([x,x],[lo,hi],color=color,lw=1.3,zorder=2)
            for y,truncated in [(lo,region['lower_grid_truncated']),(hi,region['upper_grid_truncated'])]:
                if not truncated:ax.plot(x,y,marker='_',color=color,ms=7,mew=1.2,zorder=2)
        ax.plot(x,r['rate_msun_myr'],marker=markers[r['tail_shape']],ms=6.8,
            markerfacecolor='white' if r['rate_best_on_grid_edge'] else color,markeredgecolor=color,markeredgewidth=1.3,lw=0,zorder=3)
    add_labels(ax,rows,xkey)
    if fit is not None:
        xx=np.linspace(*ax.get_xlim(),200);ax.plot(xx,fit['a']*xx+fit['b'],color='#374151',ls='--',lw=1.4,zorder=1)
        sign='+' if fit['b']>=0 else '-'
        text=rf"$y={fit['a']:.3f}x{sign}{abs(fit['b']):.3f}$"+'\n'+rf"$R^2={fit['r_squared']:.3f}$"+f"，N={fit['n']}"
        ax.text(.97,.96,text,transform=ax.transAxes,ha='right',va='top',fontsize=10,
                bbox=dict(facecolor='white',edgecolor='#D7DDE3',boxstyle='round,pad=.4',alpha=.95))

def run(out,doc,fit_linear=False):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    os.environ['MPLCONFIGDIR']=str(out/'mplcache')
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib import font_manager
    font_manager.fontManager.addfont('C:/Windows/Fonts/msjh.ttc')
    plt.rcParams.update({'font.family':'Microsoft JhengHei','axes.unicode_minus':False,'font.size':10,
                         'axes.titlesize':12,'axes.labelsize':11,'pdf.fonttype':42,'ps.fonttype':42})
    master_path=STUDY/'catalogs/master_clusters.csv';master=pd.read_csv(master_path).set_index('id')
    shapes,physical=document_tables(doc);crosschecks=[]
    for cid,values in physical.items():
        r=master.loc[cid];ok=(abs(r['mass']-values['mass'])<1 and abs(r['ecc']-values['ecc'])<=.0051 and abs(r.R_gc-values['R_gc'])<=.051)
        crosschecks.append(dict(cluster=cid,document_and_catalogue_agree_at_document_precision=bool(ok)))
        if not ok:raise ValueError('Supplied physical table differs from the frozen catalogue: '+cid)
    estimates,excluded=load_estimates();records=[]
    for cid,r in estimates.items():
        cat=master.loc[cid];profile=sorted(r['profile'],key=lambda p:p['rate_msun_myr']);x=[p['rate_msun_myr'] for p in profile];q=[p['q'] for p in profile]
        best=float(r['best']['rate_msun_myr']);regions=q_regions(x,q)
        selected=next((s for s in regions if s['lower']-1e-9<=best<=s['upper']+1e-9),None)
        if selected is None:raise ValueError('Best grid point not in diagnostic support')
        eccentricity=(cat.R_apo-cat.R_peri)/(cat.R_apo+cat.R_peri)
        if not math.isclose(cat.ecc,eccentricity,abs_tol=5.1e-6):raise ValueError('Eccentricity definition mismatch')
        records.append(dict(cluster=cid,label=LABELS.get(cid,cid),delta=0.,rate_msun_myr=best,
            q1_lower=selected['lower'],q1_upper=selected['upper'],
            diagnostic_minus=None if selected['lower_grid_truncated'] else best-selected['lower'],
            diagnostic_plus=None if selected['upper_grid_truncated'] else selected['upper']-best,
            lower_grid_truncated=selected['lower_grid_truncated'],upper_grid_truncated=selected['upper_grid_truncated'],
            rate_best_on_grid_edge=best in (min(x),max(x)),R_gc_kpc=float(cat.R_gc),mass_msun=float(cat['mass']),eccentricity=float(eccentricity),
            mass_1e5_msun=float(cat['mass']/1e5),
            r_peri_kpc=float(cat.R_peri),r_apo_kpc=float(cat.R_apo),tail_shape=shapes[cid]['shape'],
            q1_regions=regions,rate_profile=profile,sources=r['sources'],confidence_intervals_calibrated=False,
            target_posterior_computed=False,physical_rate_measurement_validated=False))
    records.sort(key=lambda r:r['R_gc_kpc'])
    fit_keys=['R_gc_kpc','mass_1e5_msun','eccentricity']
    fits={key:linear_fit(records,key) for key in fit_keys} if fit_linear else {}
    if fits:pd.DataFrame(fits.values()).to_csv(out/'linear_fits.csv',index=False,encoding='utf-8-sig')
    flat=pd.DataFrame([{k:v for k,v in r.items() if k not in ('q1_regions','rate_profile','sources')} for r in records])
    flat.to_csv(out/'plot_data_full_precision.csv',index=False,encoding='utf-8-sig')
    flat.to_csv(out/'plot_data_two_decimals.csv',index=False,encoding='utf-8-sig',float_format='%.2f')
    report=dict(created_taipei=datetime.now(timezone(timedelta(hours=8))).isoformat(),clusters=len(records),delta=0.,
        statistic='Integrated-background-and-MF likelihood; best sampled rate after profiling over sampled T',
        vertical_spans='Connected q_profile<=1 regions using piecewise-linear interpolation between computed rate nodes; not confidence or credible intervals.',
        grid_rate_spacing_msun_myr=2.,display_precision_does_not_imply_measurement_accuracy=True,
        sources={str(master_path):sha(master_path),str(Path(doc).resolve()):sha(doc)},document_crosschecks=crosschecks,
        parameter_sources=['https://people.smp.uq.edu.au/HolgerBaumgardt/globular/combined_table.txt',
            'https://people.smp.uq.edu.au/HolgerBaumgardt/globular/orbits.html'],
        eccentricity_definition='(R_apo-R_peri)/(R_apo+R_peri), from the catalogue orbit model, not recalculated in the spray MWPotential2014',
        excluded=excluded,records=records,linear_fits=fits,confidence_intervals_calibrated=False,physical_rate_measurement_validated=False)
    (out/'data_and_provenance.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    legend=[Line2D([],[],marker=m,color=c,lw=0,label=label,markersize=6) for label,c,m in [('雙側','#2878AD','o'),('S 形','#B87520','s'),('單側','#8E5196','^')]]
    legend.append(Line2D([],[],marker='o',mfc='white',mec='#53616F',lw=0,label='網格邊界最佳點',markersize=6))
    specs=[('R_gc_kpc','MLR 與銀心距離',r'銀心距離 $R_{\rm GC}$ (kpc)','mlr_vs_galactocentric_distance'),
           ('mass_1e5_msun' if fit_linear else 'mass_msun','MLR 與星團質量',
            r'星團質量 $M$ ($10^5\,M_\odot$)' if fit_linear else r'星團質量 $M$ ($M_\odot$)','mlr_vs_mass'),
           ('eccentricity','MLR 與軌道偏心率',r'軌道偏心率 $e$','mlr_vs_eccentricity')]
    note='垂直棒：'+r'$q_{\rm prof}\leq1$'+' 診斷範圍，網格間線性插值；不是 1σ 或信賴區間。\n無端帽的一端受網格截斷；形態分類依提供資料表，本次未重新驗證。'
    if fit_linear:note+='\n虛線：普通最小平方 y=ax+b；包含網格邊界點，未以診斷棒加權。'
    for key,title,xlabel,name in specs:
        fig,ax=plt.subplots(figsize=(8,6.6));fig.subplots_adjust(left=.11,right=.97,top=.82,bottom=.24)
        panel(ax,records,key,title,xlabel,fits.get(key));fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.52,.965),ncol=4,frameon=False,fontsize=9)
        fig.text(.98,.90,f'δ = 0；N = {len(records)}',ha='right',fontsize=10,color='#53616F')
        fig.text(.11,.035,note,fontsize=9,color='#53616F',va='bottom',linespacing=1.7)
        for ext in ('png','pdf','svg'):fig.savefig(out/(name+'.'+ext),dpi=220)
        plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(19,6.7));fig.subplots_adjust(left=.045,right=.99,wspace=.24,top=.82,bottom=.24)
    for ax,(key,title,xlabel,_) in zip(axes,specs):panel(ax,records,key,title,xlabel,fits.get(key))
    fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.5,.96),ncol=4,frameon=False,fontsize=10)
    fig.text(.5,.875,f'δ = 0；N = {len(records)}；同一定義的積分概似結果',ha='center',fontsize=11,color='#53616F')
    fig.text(.045,.025,note,fontsize=10,color='#53616F',va='bottom',linespacing=1.7)
    fig.savefig(out/'three_parameter_comparison.png',dpi=200);fig.savefig(out/'three_parameter_comparison.pdf');plt.close(fig)
    lines=['# δ＝0 質量流失率與星團參數比較','',f'本次納入 {len(records)} 顆星團。點為最佳已計算網格值；未計算 R/T 後驗。',
        '垂直棒是 qprof≤1 的診斷範圍，端點由相鄰網格間線性插值取得，不是 1σ、68% 或 95% 信賴區間。',
        '小數點後兩位是顯示格式，並不表示流失率測量精度為 0.01 M☉/Myr。原始流失率網格間距為 2 M☉/Myr。',
        '空心點表示最佳值落在流失率網格邊界；不畫端帽的一端表示 q≤1 範圍被計算網格截斷，不能解讀為物理上下限。',
        '若診斷範圍不連續，圖中保留分離的線段；下表列包含最佳點的區段。被網格截斷的誤差方向以「未定」標示，不以零誤差呈現。','',
        '| 星團 | MLR | 診斷 − | 診斷 ＋ | q≤1 區段 | 銀心距離 (kpc) | 質量 (M☉) | 偏心率 |',
        '|---|---:|---:|---:|---|---:|---:|---:|']
    for r in records:
        lo='未定' if r['diagnostic_minus'] is None else f"{r['diagnostic_minus']:.2f}"
        hi='未定' if r['diagnostic_plus'] is None else f"{r['diagnostic_plus']:.2f}"
        lines.append(f"| {r['label']} | {r['rate_msun_myr']:.2f} | {lo} | {hi} | {r['q1_lower']:.2f}–{r['q1_upper']:.2f} | {r['R_gc_kpc']:.2f} | {r['mass_msun']:,.0f} | {r['eccentricity']:.2f} |")
    lines+=['','## 參數來源與解讀','',
        '採本機保存的 Baumgardt 目錄版本，與本專案設定一致，沒有改用即時網站更新值。所提供《Claude資料表.docx》的 29 顆質量、偏心率及銀心距離已依文件顯示精度核對。',
        '尾巴形態顏色依文件表5，僅作資料分組，本次未重新驗證潮汐尾偵測或形態分類。',
        '銀心距離為目前的三維銀心距離，不是近銀心距離。偏心率採 (Rapo−Rperi)/(Rapo+Rperi)，源自目錄採用的軌道模型。',
        '圖中未畫橫軸誤差，不代表質量、距離與軌道沒有不確定性。這些參數亦參與模型設定，圖上關聯不能直接當作獨立因果證據。',
        ('依要求對三個橫軸各擬合 y=ax+b，使用普通最小平方且不以 q 診斷棒加權。質量以 10^5 M☉ 為單位，未對質量取對數；直線亦包含網格邊界最佳點，僅作描述性擬合，不宣稱物理定律、因果關係或顯著性。'
         if fit_linear else '此次不擬合趨勢線，也不將診斷棒用作已知標準差進行加權迴歸。'),
        'NGC362 缺少同一定義的 comparison，未將最高單一情境的 scores 結果混入。其餘尚未完成 δ=0 比較的星團也未納入。','',
        '來源：[Baumgardt 質量與結構目錄](https://people.smp.uq.edu.au/HolgerBaumgardt/globular/combined_table.txt)、[軌道目錄與模型說明](https://people.smp.uq.edu.au/HolgerBaumgardt/globular/orbits.html)。',
        '統計區間概念可參考 [PDG 統計綜述](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-statistics.pdf)；本圖未宣稱具有頻率學涵蓋率。']
    if fits:
        lines+=['','## y=ax+b 描述性擬合','','| x | a | b | R² | N |','|---|---:|---:|---:|---:|']
        for key,fit in fits.items():lines.append(f"| {key} | {fit['a']:.6f} | {fit['b']:.6f} | {fit['r_squared']:.6f} | {fit['n']} |")
        lines+=['','銀心距離 x 單位為 kpc；質量 x 為 M/(10^5 M☉)；偏心率無因次。若質量改用 M☉，a 須再除以 10^5，b 與 R² 不變。','係數以獨立的中心化 OLS 公式交叉核對；未計算係數信賴區間或相關性 p 值。']
    (out/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(dict(clusters=len(records),excluded=excluded,out=str(out),disconnected_q_regions=[r['cluster'] for r in records if len(r['q1_regions'])>1]),ensure_ascii=True))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--doc',type=Path,required=True)
    p.add_argument('--fit-linear',action='store_true');a=p.parse_args();run(a.out,a.doc,a.fit_linear)
