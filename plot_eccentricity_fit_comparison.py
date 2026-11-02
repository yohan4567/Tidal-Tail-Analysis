"""Compare linear and a/(1-b*e) fits to an audited, fixed MLR sample."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares, minimize_scalar

from plot_selected_mlr import BLUE,ols,panel,sha


def rational_fit(rows):
    x=np.array([r['eccentricity'] for r in rows],float)
    y=np.array([r['rate_msun_myr'] for r in rows],float)
    if not (np.isfinite(x).all() and np.isfinite(y).all() and (x>0).all() and (x<1).all() and (y>0).all()):
        raise ValueError('Need positive finite MLR and 0 < e < 1')
    if np.unique(x).size<2:raise ValueError('Need distinct eccentricities')
    upper=1/float(x.max())
    def profile(b):
        g=1/(1-b*x);a=float(g@y/(g@g));return float(np.sum((y-a*g)**2))
    # Profile out a analytically, and cross-check with joint optimization.
    p=minimize_scalar(profile,bounds=(-10000.,upper-1e-9),method='bounded',options={'xatol':1e-12})
    b=float(p.x);g=1/(1-b*x);a=float(g@y/(g@g));sse=profile(b)
    checks=[least_squares(lambda v:v[0]/(1-v[1]*x)-y,[float(y.mean()),start],
        bounds=([0.,-np.inf],[np.inf,upper-1e-9]),max_nfev=10000,
        ftol=1e-12,xtol=1e-12,gtol=1e-12) for start in (-100.,-10.,-1.,0.,.5)]
    independent=min(float(v.fun@v.fun) for v in checks)
    if not p.success or not np.isclose(sse,independent,rtol=1e-9,atol=1e-9):
        raise ValueError('Profile and joint rational optimizations disagree')
    # Check that the optimum is finite rather than the b -> -infinity limit.
    inv=1/x;c=float(inv@y/(inv@inv));limit_sse=float(np.sum((y-c*inv)**2))
    if not sse<limit_sse-1e-8 or b<=-9999:raise ValueError('No supported finite negative-b optimum')
    if (1-b*x<=0).any():raise ValueError('Denominator is not positive over the sample')
    return dict(model='a_over_one_minus_b_e',equation='y = a / (1 - b * e)',a=a,b=b,
        r_squared=1-sse/float(np.sum((y-y.mean())**2)),residual_sum_squares=sse,
        n=len(rows),fit_clusters=[r['cluster'] for r in rows],parameters=2,
        domain='a > 0, b < 1/max(observed eccentricity); b may be negative',
        b_negative_infinity_limit_SSE=limit_sse,
        minimum_sample_denominator=float(np.min(1-b*x)),
        independent_joint_fit_check_passed=True,confidence_intervals_calibrated=False)


def run(source,out):
    source,out=Path(source).resolve(),Path(out).resolve()
    if out.exists():raise FileExistsError('Use a fresh output version')
    data=json.loads(source.read_text(encoding='utf-8'));rows=data['records']
    if len(rows)<3 or len({r['cluster'] for r in rows})!=len(rows):raise ValueError('Invalid sample membership')
    if not all(r['delta']==0 and r['method']=='background_and_MF_integrated' for r in rows):
        raise ValueError('Do not mix offsets or estimators')
    fits={'linear':ols(rows,'eccentricity'),'rational':rational_fit(rows)}
    fits['linear'].update(model='linear',equation='y = a * e + b',parameters=2)
    x=np.array([r['eccentricity'] for r in rows]);y=np.array([r['rate_msun_myr'] for r in rows])
    for name,f in fits.items():
        prediction=f['a']*x+f['b'] if name=='linear' else f['a']/(1-f['b']*x)
        sse=float(np.sum((y-prediction)**2));sst=float(np.sum((y-y.mean())**2))
        if not np.isclose(f['r_squared'],1-sse/sst):raise ValueError('R squared check failed')
        f.update(residual_sum_squares=sse,total_sum_squares=sst,
                 weighting='unweighted least squares in original MLR units')
    out.mkdir(parents=True);(out/'data').mkdir();os.environ['MPLCONFIGDIR']=str(out/'data/mplcache')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.lines import Line2D
    font_manager.fontManager.addfont('C:/Windows/Fonts/msjh.ttc')
    plt.rcParams.update({'font.family':['Microsoft JhengHei','DejaVu Sans'],'axes.unicode_minus':False,
        'font.size':10,'axes.titlesize':13,'axes.labelsize':12})
    colors={'linear':'#59616C','rational':'#B26927'};styles={'linear':'--','rational':'-'}
    exported=[]
    for names,stem in [(('linear','rational'),'mlr_vs_eccentricity_two_fits'),
                       (('linear',),'mlr_vs_eccentricity_linear'),(('rational',),'mlr_vs_eccentricity_rational')]:
        fig,ax=plt.subplots(figsize=(9.4,8.2));fig.subplots_adjust(left=.12,right=.98,bottom=.23,top=.82)
        # Keep the existing data-point, label and q-diagnostic rendering.
        panel(ax,rows,'eccentricity','MLR 與軌道偏心率：擬合比較',r'軌道偏心率 $e$',fits['linear'])
        ax.lines[-1].remove();ax.texts[-1].remove()
        xx=np.linspace(float(x.min()),float(x.max()),400)
        handles=[Line2D([],[],marker='o',color=BLUE,lw=0,label=f'{len(rows)} 顆星團'),
            Line2D([],[],marker='o',mfc='white',mec='#59616C',lw=0,label='網格邊界最佳值；無端帽＝截斷')]
        labels=[]
        for name in names:
            f=fits[name]
            if name=='linear':
                predicted=f['a']*xx+f['b'];sign='+' if f['b']>=0 else '-'
                equation=rf"$y={f['a']:.3f}e{sign}{abs(f['b']):.3f}$";title='直線'
            else:
                predicted=f['a']/(1-f['b']*xx);sign='-' if f['b']>=0 else '+'
                equation=rf"$y=\frac{{{f['a']:.3f}}}{{1{sign}{abs(f['b']):.3f}e}}$";title='a/(1-be)'
            ax.plot(xx,predicted,color=colors[name],ls=styles[name],lw=1.6,zorder=1)
            handles.append(Line2D([],[],color=colors[name],ls=styles[name],lw=1.6,label=title))
            labels.append(equation+rf"，$R^2={f['r_squared']:.3f}$")
        ax.text(.97,.97,'\n'.join(labels),transform=ax.transAxes,ha='right',va='top',fontsize=11,
            linespacing=2.1,bbox=dict(facecolor='white',edgecolor='#D7DDE3',boxstyle='round,pad=.5',alpha=.96))
        fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.55,.985),ncol=2,frameon=False,fontsize=10)
        fig.text(.98,.87,f'offset＝0；繪圖及擬合 N＝{len(rows)}',ha='right',fontsize=10,color='#53616F')
        fig.text(.12,.04,r'診斷線段：$q_{\rm prof}\leq1$，不是標準差或信賴區間。'
            '\n兩模型均在原始 MLR 尺度作未加權最小平方，未將診斷棒當作標準差加權。'
            '\na/(1-be) 的 b 可為負值；分母在資料範圍內保持正值，不允許穿越奇點。',
            fontsize=9,color='#53616F',linespacing=1.7,va='bottom')
        for ext in ('png','pdf'):
            path=out/(stem+'.'+ext);fig.savefig(path,dpi=240,facecolor='white');exported.append(path.name)
        plt.close(fig)
    (out/'data/provenance.json').write_text(json.dumps(dict(source=str(source),source_sha256=sha(source),
        code_sha256=sha(Path(__file__)),fits=fits,records=rows,exports={p:sha(out/p) for p in exported},
        target_posterior_used=False,confidence_intervals_calibrated=False,
        interpretation='Descriptive fits to fixed conditional sampled estimates; model comparison is not a physical validation.'),ensure_ascii=False,indent=2),encoding='utf-8')
    pd.DataFrame(fits.values()).to_csv(out/'data/fit_parameters.csv',index=False,encoding='utf-8-sig')
    print(json.dumps(dict(out=str(out),n=len(rows),fits=fits),ensure_ascii=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=Path('outputs/mlr_three_comparisons_14_20261102_v1/data/provenance.json'))
    p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.source,a.out)
