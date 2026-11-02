"""Plot existing calibration counts, without changing selections or thresholds."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from provenance import sha256_file
from selection_model import normalize_gaia, science_available


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def project(ra,dec,ra0,dec0):
    a,d,a0,d0=map(np.deg2rad,(ra,dec,ra0,dec0))
    denominator=np.sin(d0)*np.sin(d)+np.cos(d0)*np.cos(d)*np.cos(a-a0)
    return (np.rad2deg(np.cos(d)*np.sin(a-a0)/denominator),
            np.rad2deg((np.cos(d0)*np.sin(d)-np.sin(d0)*np.cos(d)*np.cos(a-a0))/denominator))


def run(bundle,selection,out):
    bundle,selection,out=map(lambda p:Path(p).resolve(),(bundle,selection,out))
    if out.exists():raise FileExistsError('Use a fresh diagnostic output directory')
    meta=read(bundle/'quality/quality.json');report=read(selection/'selection.json')
    parent_path=bundle/'parent_spatial_calibration.csv.gz';npz=bundle/'quality/quality.npz'
    if sha256_file(npz)!=meta['table_sha256']:raise ValueError('Stored count table checksum changed')
    if sha256_file(parent_path)!=meta['parent_sha256']:raise ValueError('Parent catalogue checksum changed')
    if meta['selection_sha256']!=report['selection_file_sha256']:raise ValueError('Selection identities differ')
    if [meta['center_ra_deg'],meta['center_dec_deg']]!=[report['center_ra_deg'],report['center_dec_deg']]:
        raise ValueError('Centres differ')
    with np.load(npz) as z:a={k:z[k] for k in z.files}
    d=normalize_gaia(pd.read_csv(parent_path,low_memory=False))
    x,y=project(d.ra.to_numpy(),d.dec.to_numpy(),meta['center_ra_deg'],meta['center_dec_deg'])
    axes=[a[k] for k in ('x_edges','y_edges','g_edges','c_edges')]
    xyz=np.column_stack((x,y,d.phot_g_mean_mag));available=science_available(d)
    xyzc=np.column_stack((xyz,d.g_rp))
    # Reconstruct integer counts from the actual catalogue before plotting.
    for key,features,bins in [('n_parent',xyz,axes[:3]),('k_available',xyz[available],axes[:3]),
                              ('n_available',xyzc[available],axes)]:
        if not np.array_equal(np.histogramdd(features,bins=bins)[0],a[key]):
            raise ValueError('Catalogue count reconstruction differs: '+key)
    minimum=meta['min_parent'];counts=a['n_available'];bad=counts<minimum
    rows=[]
    for ix,iy,ig,ic in np.ndindex(counts.shape):
        rows.append(dict(x_bin=ix,y_bin=iy,g_bin=ig,colour_bin=ic,
            x_min=axes[0][ix],x_max=axes[0][ix+1],y_min=axes[1][iy],y_max=axes[1][iy+1],
            g_min=axes[2][ig],g_max=axes[2][ig+1],g_rp_min=axes[3][ic],g_rp_max=axes[3][ic+1],
            available_reference_stars=int(counts[ix,iy,ig,ic]),minimum=minimum,
            below_minimum=bool(bad[ix,iy,ig,ic]),quality_passing_stars=int(a['k_quality'][ix,iy,ig,ic])))
    out.mkdir(parents=True);os.environ['MPLCONFIGDIR']=str(out/'mplcache')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.patches import Circle,Rectangle,Patch
    font_manager.fontManager.addfont('C:/Windows/Fonts/msjh.ttc')
    plt.rcParams.update({'font.family':['Microsoft JhengHei','DejaVu Sans'],'axes.unicode_minus':False,
        'font.size':10,'pdf.fonttype':42,'ps.fonttype':42})
    radius=np.rad2deg(np.tan(np.deg2rad(meta['parent_metadata']['radius_deg'])))
    core=meta['parent_metadata']['client_spatial_exclusions']['core_radius_deg']
    cluster=report['cluster'];red='#B44439';blue='#2C6484'

    def sky(ax,title,small=False):
        ax.set_facecolor('#F3F4F6')
        disk=Circle((0,0),radius,facecolor='white',edgecolor='#64748B',lw=.9,zorder=0)
        ax.add_patch(disk)
        ax.axvline(0,color='#64748B',lw=.7,zorder=1);ax.axhline(0,color='#64748B',lw=.7,zorder=1)
        ax.set_xlim(radius*1.04,-radius*1.04);ax.set_ylim(-radius*1.04,radius*1.04);ax.set_aspect('equal')
        ax.set_title(title,fontsize=10 if small else 13,pad=10)
        ax.set_xticks([-2,0,2]);ax.set_yticks([-2,0,2]);ax.tick_params(labelsize=8 if small else 10)
        if not small:ax.set_xlabel('切平面東向偏移 ξ（度；東在左）');ax.set_ylabel('切平面北向偏移 η（度）')
        for spine in ax.spines.values():spine.set_visible(False)
        return disk

    def masks(ax):
        ax.add_patch(Circle((0,0),core,facecolor='#D8DCE1',edgecolor='#64748B',lw=.7,zorder=4))
        ax.plot(0,0,'+',ms=6,color='#334155',zorder=5)

    def cellmap(ax,n,disk,highlight=True,small=False):
        for ix,iy in np.ndindex(n.shape):
            lo,hi=axes[0][ix:ix+2];bottom,top=axes[1][iy:iy+2]
            isbad=highlight and n[ix,iy]<minimum
            color=red if isbad else blue
            patch=Rectangle((lo,bottom),hi-lo,top-bottom,facecolor='#FBE7E4' if isbad else '#EDF4F8',edgecolor='none',zorder=.5)
            patch.set_clip_path(disk);ax.add_patch(patch)
            label=f'{int(n[ix,iy]):,}'+('\n未達標' if isbad else '')
            ax.text((lo+hi)/2,(bottom+top)/2,label,ha='center',va='center',
                fontsize=11 if small else 22,color=color,fontweight='bold',zorder=6)
        masks(ax)

    def save(fig,name):
        for ext in ('png','pdf','svg'):fig.savefig(out/(name+'.'+ext),dpi=220,facecolor='white')
        plt.close(fig)

    # Overview: total counts alone conceal sparse magnitude/colour strata.
    fig,axs=plt.subplots(1,2,figsize=(12.8,7.5))
    fig.subplots_adjust(left=.07,right=.97,bottom=.24,top=.78,wspace=.28)
    disk=sky(axs[0],'四個校正天區：母樣本總星數')
    axs[0].scatter(x,y,s=.12,alpha=.18,c='#657487',rasterized=True,lw=0,zorder=.8)
    totals=a['n_parent'].sum(axis=2);fail_per_sky=bad.sum(axis=(2,3))
    for ix,iy in np.ndindex(totals.shape):
        cx=axes[0][ix:ix+2].mean();cy=axes[1][iy:iy+2].mean()
        axs[0].text(cx,cy,f'N＝{int(totals[ix,iy]):,}\n{int(fail_per_sky[ix,iy])}／{counts.shape[2]*counts.shape[3]} 格未達標',
            ha='center',va='center',fontsize=13,color='#23364A',zorder=6,
            bbox=dict(facecolor='white',alpha=.94,edgecolor='#CDD5DD',boxstyle='round,pad=.5'))
    masks(axs[0])
    # Select the stratum with the most unsupported sky cells; ties use least support.
    ig,ic=min(np.ndindex(counts.shape[2:]),key=lambda j:(-int(bad[:,:,j[0],j[1]].sum()),int(counts[:,:,j[0],j[1]].sum())))
    g0,g1=axes[2][ig:ig+2];c0,c1=axes[3][ic:ic+2]
    example=dict(g_range=[float(g0),float(g1)],g_rp_range=[float(c0),float(c1)],counts_xy=counts[:,:,ig,ic].tolist())
    disk=sky(axs[1],f'稀疏分層：{g0:g} ≤ G < {g1:g}\n{c0:g} ≤ G-RP < {c1:g}，可用參考星數')
    cellmap(axs[1],counts[:,:,ig,ic],disk)
    fig.suptitle(f'{cluster} 天區與校正樣本支持度',fontsize=17,y=.96)
    fig.text(.5,.89,f'品質校正門檻：每個天區 × 星等 × 顏色格至少 {minimum} 顆可用參考星；{int(bad.sum())}／{bad.size} 格未達標',ha='center',fontsize=11)
    fig.text(.07,.045,'紅色：可用參考星數低於門檻；灰色中央圓：校正母樣本排除的擁擠核心。\n'
        '圖中為 Gaia DR3 校正場星，包含背景星，並非星團成員數。G < 12 採全量，12 ≤ G < 21 採 10% 獨立樣本。\n'
        '左圖總數不能取代分層檢查；右圖僅展示最稀疏的分層，全部分層另見完整圖。',fontsize=10,color='#53616F',linespacing=1.7)
    save(fig,'NGC6535_sky_support_overview')

    fig,axs=plt.subplots(counts.shape[3],counts.shape[2],figsize=(16,7.7),squeeze=False)
    fig.subplots_adjust(left=.055,right=.98,bottom=.18,top=.78,wspace=.16,hspace=.3)
    for ic in range(counts.shape[3]):
        for ig in range(counts.shape[2]):
            ax=axs[ic,ig];g0,g1=axes[2][ig:ig+2];c0,c1=axes[3][ic:ic+2]
            disk=sky(ax,f'{g0:g} ≤ G < {g1:g}\n{c0:g} ≤ G-RP < {c1:g}',small=True)
            cellmap(ax,counts[:,:,ig,ic],disk,small=True)
    fig.suptitle(f'{cluster} 全部品質校正格：可用參考星數',fontsize=17,y=.96)
    fig.legend(handles=[Patch(facecolor='#FBE7E4',edgecolor=red,label=f'未達標：N < {minimum}'),
                        Patch(facecolor='#EDF4F8',edgecolor=blue,label=f'星數達標：N ≥ {minimum}')],
        loc='upper center',bbox_to_anchor=(.5,.915),ncol=2,frameon=False)
    fig.text(.055,.035,f'每小圖為相同四天區；橫向為東向偏移 ξ（東在左），縱向為北向偏移 η，單位為度。\n'
        f'共 {int(bad.sum())}／{bad.size} 個品質校正格低於 {minimum} 顆；星數達標不表示所有其他科學驗證已通過。\n'
        '數字為品質篩選前、必要欄位可用的校正場星數；灰色中央圓為校正母樣本的核心排除區。',fontsize=10,color='#53616F',linespacing=1.6)
    save(fig,'NGC6535_all_quality_cells')
    pd.DataFrame(rows).to_csv(out/'quality_cell_counts.csv',index=False,encoding='utf-8-sig')
    summary=dict(cluster=cluster,min_reference_stars=minimum,total_quality_cells=int(bad.size),
        unsupported_quality_cells=int(bad.sum()),availability_cells=int(a['n_parent'].size),
        unsupported_availability_cells=int((a['n_parent']<minimum).sum()),
        retained_parent_rows=len(d),parent_rows_in_grid=int(a['n_parent'].sum()),
        total_parent_counts_xy=totals.tolist(),unsupported_quality_cells_xy=fail_per_sky.tolist(),
        example_stratum=example,
        count_reconstruction_passed=True,
        sources={str(p):sha256_file(p) for p in (npz,bundle/'quality/quality.json',parent_path,selection/'selection.json')},
        interpretation='Calibration field-star support only; not globular-cluster member counts, tail significance or a changed calibration.')
    (out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'README.md').write_text(f'# {cluster} 校正天區星數\n\n'
        f'品質校正 {bad.size} 格中有 {int(bad.sum())} 格的可用參考星數低於 {minimum}。'
        '四個空間格與原始 G、G-RP 邊界皆未更動。數字為校正場星，不是星團成員。\n\n'
        'overview 為天區總數與稀疏分層示例；all_quality_cells 為全部分層。'
        '核心依母樣本原規則排除；母樣本超出矩形校正網格的少量星不計入格內總數。'
        '原資料、校正、門檻及正在執行的分析均未修改。\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k!='sources'},ensure_ascii=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bundle',type=Path,required=True);p.add_argument('--selection',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.bundle,a.selection,a.out)
