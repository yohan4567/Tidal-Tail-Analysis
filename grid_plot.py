"""Plot sampled T/rate cells without interpolating missing grid nodes."""
import os
from pathlib import Path
import numpy as np
import pandas as pd


def edges(values):
    v=np.asarray(values,float)
    if len(v)==1:return np.array([v[0]-.5,v[0]+.5])
    mid=(v[:-1]+v[1:])/2
    return np.r_[v[0]-(mid[0]-v[0]),mid,v[-1]+(v[-1]-mid[-1])]


def plot_grid(rows,path,confidence=None):
    os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).parent/'outputs/.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    f=pd.DataFrame(rows).groupby(['duration_myr','rate_msun_myr'],as_index=False).q.min()
    times=sorted(f.duration_myr.unique());rates=sorted(f.rate_msun_myr.unique())
    z=np.full((len(rates),len(times)),np.nan)
    for row in f.itertuples():z[rates.index(row.rate_msun_myr),times.index(row.duration_myr)]=row.q
    fig,ax=plt.subplots(figsize=(9,6),layout='constrained')
    mesh=ax.pcolormesh(edges(times)/1000,edges(rates),np.ma.masked_invalid(z),cmap='viridis_r',shading='flat')
    fig.colorbar(mesh,ax=ax,label='q = 2 (NLL - minimum NLL); smaller is better')
    ax.set_xticks(np.asarray(times)/1000); ax.set_yticks(rates)
    ax.scatter(f.duration_myr/1000,f.rate_msun_myr,s=12,c='black',alpha=.65,label='Computed grid nodes')
    if len(f)<=60:
        for row in f.itertuples():
            ax.annotate(f'{row.q:.2f}',(row.duration_myr/1000,row.rate_msun_myr),xytext=(5,6),textcoords='offset points',fontsize=8,color='black',bbox=dict(facecolor='white',alpha=.7,edgecolor='none',pad=1))
    best=f.loc[f.q.idxmin()];ax.scatter([best.duration_myr/1000],[best.rate_msun_myr],marker='*',s=150,c='red',edgecolors='white',label='Best sampled node')
    if confidence is not None:
        accepted=[r for r in confidence if r['accepted95']]
        if accepted:ax.scatter([r['duration_myr']/1000 for r in accepted],[r['rate_msun_myr'] for r in accepted],s=120,facecolors='none',edgecolors='white',linewidths=1.5,label='In conditional 95% joint set')
    ax.set(xlabel='Lookback window T (Gyr)',ylabel='Escape rate (solar masses/Myr)',
           title='Conditional T–MLR grid; population profiled\nDiscrete nodes; physical validation pending')
    ax.legend(fontsize=8);fig.savefig(path,dpi=170);plt.close(fig)
