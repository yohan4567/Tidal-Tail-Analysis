"""Population-conditioned rate diagnostics without changing any fitted likelihood.

Post-processing only. Existing grid/response source identities are not changed.
"""
import argparse
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import time
import numpy as np
import pandas as pd
from scipy.special import xlogy
from provenance import sha256_file
from joint_grid import aggregate_pixels


def population_summary(rows):
    groups={}
    for row in rows:groups.setdefault(row.get('scenario','unspecified'),[]).append(row)
    if not groups:raise ValueError('No grid candidates')
    global_min=min(r['negative_log_likelihood'] for r in rows)
    result=[]
    for scenario,items in sorted(groups.items()):
        best=min(items,key=lambda r:r['negative_log_likelihood']);rates=sorted({r['rate_msun_myr'] for r in items})
        times=sorted({r['duration_myr'] for r in items})
        rate=float(best['rate_msun_myr']);duration=float(best['duration_myr'])
        mass=rate*duration
        result.append(dict(scenario=scenario,rate_msun_myr=rate,duration_myr=duration,
            q_relative_to_global_min=2*(best['negative_log_likelihood']-global_min),
            expected_signal_stars=best.get('expected_signal_stars'),
            represented_escaped_mass_in_window_msun=mass,
            selected_stars_per_escaped_msun=None if best.get('expected_signal_stars') is None else best['expected_signal_stars']/mass,
            rate_at_grid_boundary=rate in (rates[0],rates[-1]),time_at_grid_boundary=duration in (times[0],times[-1])))
    rates=[r['rate_msun_myr'] for r in result]
    return dict(primary_rate_msun_myr=None,physical_rate_measurement_validated=False,
        population_normalization_status='not_externally_calibrated',conditional_scenario_results=result,
        conditional_best_rate_range_msun_myr=[min(rates),max(rates)],conditional_best_rate_spread_factor=max(rates)/min(rates),
        range_is_confidence_interval=False,
        interpretation='Scenario sensitivity only; the best-rate span is neither a confidence interval nor a calibrated systematic-error bound. Boundary solutions can extend beyond it.',
        visibility_definition='Selected stars / (rate * lookback duration). Includes dynamics, footprint, time window, photometry and selection; not a universal population-only factor.',
        required_external_information='Escaping-population luminosity/mass-function and unseen-mass constraints with an explicit selection model; a core mass function alone is not an escaping-population measurement.')


def _label(scenario):
    try:
        s=json.loads(scenario)
        return ', '.join(f'{k}={v:g}' if isinstance(v,(int,float)) else f'{k}={v}' for k,v in (s or {}).items()) or 'unspecified'
    except (ValueError,TypeError):return str(scenario)


def selected_arrays(score,row):
    key=lambda r:(r['duration_myr'],r['rate_msun_myr'],r['scenario'])
    diag=next(d for d in score['diagnostics'] if key(d)==key(row))
    states=[];n=b=description=None
    for item in diag['inputs']:
        path=Path(item['path'])
        if sha256_file(path)!=item['sha256']:raise ValueError('Scored template changed: '+str(path))
        with np.load(path,allow_pickle=False) as data:
            nn=data['counts'].ravel();bb=data['background'].ravel();ss=data['signal_draws']
            desc=json.loads(str(data['selection_description'].item()))
        if n is None:n,b,description=nn,bb,desc
        elif not np.array_equal(n,nn) or not np.array_equal(b,bb):raise ValueError('Candidate inputs do not share observations/background')
        states.append(ss)
    signal=np.mean(states,axis=0)
    index=score.get('pixel_groups')
    if index is not None:
        n=aggregate_pixels(n,index);b=aggregate_pixels(b,index);signal=aggregate_pixels(signal,index)
    mean_signal=signal.mean(axis=0)
    if not np.isclose(mean_signal.sum(),row['expected_signal_stars'],rtol=1e-10,atol=1e-9):raise ValueError('Scored signal reconstruction mismatch')
    beta=float(row['background_scale']);expected=mean_signal+beta*b
    residual=np.divide(n-expected,np.sqrt(expected),out=np.zeros_like(expected),where=expected>0)
    impossible=(expected==0)&(n>0);residual[impossible]=np.inf
    deviance=2*(xlogy(n,np.divide(n,expected,out=np.ones_like(expected),where=expected>0))-(n-expected))
    deviance[impossible]=np.inf
    pearson=np.divide((n-expected)**2,expected,out=np.zeros_like(expected),where=expected>0);pearson[impossible]=np.inf
    return dict(observed=n,background=beta*b,signal=mean_signal,expected=expected,pearson_residual=residual,
        mean_prediction_poisson_deviance=float(deviance.sum()),mean_prediction_pearson_sum=float(pearson.sum()),
        raw_background=b,description=description,inputs=diag['inputs'])


def _coordinates(selection,index):
    with np.load(Path(selection)/'maps.npz',allow_pickle=False) as data:
        ix,iy=np.nonzero(data['fitmask']);edges=data['edges'];centres=(edges[:-1]+edges[1:])/2
    x,y=centres[ix],centres[iy]
    if index is not None:
        if len(index)!=len(x):raise ValueError('Pixel grouping/selection-map mismatch')
        _,inv=np.unique(index,return_inverse=True);count=np.bincount(inv)
        x=np.bincount(inv,weights=x)/count;y=np.bincount(inv,weights=y)/count
    return x,y


def run(grid,out=None,selection=None,compare_rates=(5.,31.)):
    grid=Path(grid).resolve();score_path=grid/'scores/joint_grid.json'
    score=json.loads(score_path.read_text(encoding='utf-8'));summary=population_summary(score['candidates'])
    out=Path(out).resolve() if out else grid/'population_diagnostics_v1'
    out.mkdir(parents=True,exist_ok=False)  # Atomic reservation: production and watcher must not both write.
    rows=pd.DataFrame(score['candidates']);selected=[]
    for _,group in rows.groupby('scenario'):selected.append(group.loc[group.negative_log_likelihood.idxmin()].to_dict())
    for rate in compare_rates:
        group=rows[np.isclose(rows.rate_msun_myr,rate,rtol=0,atol=1e-10)]
        if len(group):
            candidate=group.loc[group.negative_log_likelihood.idxmin()].to_dict()
            if not any(all(candidate[k]==r[k] for k in ('scenario','rate_msun_myr','duration_myr')) for r in selected):selected.append(candidate)
    reconstructed=[selected_arrays(score,r) for r in selected]
    if selection:
        report=json.loads((Path(selection)/'selection.json').read_text(encoding='utf-8'))
        if not reconstructed or reconstructed[0]['description']['observation_identity']!=score['observation_identity']:
            raise ValueError('Observation identity mismatch')
        # The template response provenance contains the exact frozen selection hash.
        for item in reconstructed[0]['inputs']:
            meta=json.loads(Path(item['path']).with_name('response.json').read_text(encoding='utf-8'))
            path=(Path(selection)/'selection.json').resolve()
            if meta['input_sha256'].get(str(path))!=sha256_file(path):raise ValueError('Selection provenance mismatch')
        x,y=_coordinates(selection,score.get('pixel_groups'))
    else:x=np.arange(len(reconstructed[0]['observed']));y=np.zeros_like(x)
    frame=pd.DataFrame(summary['conditional_scenario_results']);frame.to_csv(out/'scenario_best.csv',index=False)
    profiles=[]
    for scenario,g in rows.groupby('scenario'):
        reference=g.negative_log_likelihood.min()
        for rate,sub in g.groupby('rate_msun_myr'):
            best=sub.loc[sub.negative_log_likelihood.idxmin()]
            profiles.append(dict(scenario=scenario,rate_msun_myr=rate,duration_myr=best.duration_myr,
                q_within_scenario=2*(best.negative_log_likelihood-reference),q_relative_to_global_min=best.q))
    pd.DataFrame(profiles).to_csv(out/'conditional_rate_profiles.csv',index=False)
    visibility=rows.copy();visibility['represented_escaped_mass_in_window_msun']=visibility.rate_msun_myr*visibility.duration_myr
    visibility['selected_stars_per_escaped_msun']=visibility.expected_signal_stars/visibility.represented_escaped_mass_in_window_msun
    visibility.to_csv(out/'population_visibility_grid.csv',index=False)
    comparisons=[]
    for i,(row,data) in enumerate(zip(selected,reconstructed)):
        name=f'candidate_{i}'
        np.savez_compressed(out/(name+'.npz'),x=x,y=y,**{k:data[k] for k in ('observed','background','signal','expected','pearson_residual')})
        comparisons.append(dict(**row,artifact=name,mean_prediction_poisson_deviance=data['mean_prediction_poisson_deviance'],
            mean_prediction_pearson_sum=data['mean_prediction_pearson_sum'],
            inputs=data['inputs'],residual_interpretation='Descriptive residuals of the calibration-mean prediction, not a significance map or the fitted mixture-likelihood statistic.'))
    summary.update(schema_version=1,score_sha256=sha256_file(score_path),source_sha256=sha256_file(Path(__file__)),
        created_utc=datetime.now(timezone.utc).isoformat(),candidate_comparisons=comparisons,
        changed_likelihood_weights=False,external_population_constraint_applied=False,
        residuals_are_calibrated_significance=False)
    os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).parent/'outputs/.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
    for scenario,g in pd.DataFrame(profiles).groupby('scenario'):
        ax.plot(g.rate_msun_myr,g.q_within_scenario,'o-',label=_label(scenario))
    ax.set(xlabel='Escape rate (solar masses/Myr)',ylabel='q relative to each scenario minimum',
        title='Population-conditioned profiles; no calibrated confidence thresholds')
    ax.legend(fontsize=8);fig.savefig(out/'population_profiles.png',dpi=150);plt.close(fig)
    times=sorted(visibility.duration_myr.unique());reference_window=min(times,key=lambda t:abs(t-max(times)/2))
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for scenario,g in visibility[visibility.duration_myr==reference_window].groupby('scenario'):
        g=g.sort_values('rate_msun_myr');label=_label(scenario)
        axes[0].plot(g.rate_msun_myr,g.expected_signal_stars,'o-',label=label)
        axes[1].plot(g.rate_msun_myr,g.selected_stars_per_escaped_msun,'o-',label=label)
    axes[0].set(xlabel='Escape rate (solar masses/Myr)',ylabel='Expected selected stars',title=f'Common window: {reference_window:g} Myr')
    axes[1].set(xlabel='Escape rate (solar masses/Myr)',ylabel='Selected stars / escaped solar mass',title='Visibility includes dynamics, field and selection')
    axes[0].legend(fontsize=7);fig.savefig(out/'population_visibility.png',dpi=150);plt.close(fig)
    count_vmax=max(float(d[k].max()) for d in reconstructed for k in ('observed','background','expected'))
    signal_vmax=max(float(d['signal'].max()) for d in reconstructed)
    for row,data,entry in zip(selected,reconstructed,comparisons):
        fig,axes=plt.subplots(2,3,figsize=(13,8),layout='constrained')
        for ax,key,title in zip(axes.flat,('observed','background','signal','expected','pearson_residual'),
            ('Observed','Fitted background','Expected signal','Total expected','Descriptive Pearson residual')):
            value=data[key];kw=dict(cmap='coolwarm',vmin=-4,vmax=4) if key=='pearson_residual' else dict(cmap='viridis',vmin=0,vmax=(signal_vmax if key=='signal' else count_vmax))
            im=ax.scatter(x,y,c=value,marker='s',s=14,**kw);fig.colorbar(im,ax=ax,label='Descriptive standardized residual' if key=='pearson_residual' else 'Stars / fitted pixel')
            ax.set(title=title,xlabel='Gnomonic east (deg)' if selection else 'Fitted pixel index',ylabel='Gnomonic north (deg)' if selection else '')
            ax.set_aspect('equal' if selection else 'auto')
        axes.flat[-1].axis('off');axes.flat[-1].text(0,1,f"Rate: {row['rate_msun_myr']:g} solar masses/Myr\nT: {row['duration_myr']:g} Myr\nq(global): {row['q']:.3f}\n"+_label(row['scenario']).replace(', ','\n')+'\n\nNo calibrated p-value or confidence interval.\nResiduals use mean calibration prediction.',va='top',wrap=True,fontsize=9)
        fig.savefig(out/(entry['artifact']+'_residuals.png'),dpi=140);plt.close(fig)
    lines=['# Population-conditioned escape-rate results','',
        '**No validated absolute MLR is reported.** The escaping population normalization has not been externally calibrated.','',
        '| Fixed scenario | Best sampled MLR (Msun/Myr) | T (Myr) | q relative to global best | Selected stars | Rate boundary |',
        '|---|---:|---:|---:|---:|---|']
    for r in summary['conditional_scenario_results']:
        lines.append(f"| {_label(r['scenario'])} | {r['rate_msun_myr']:g} | {r['duration_myr']:g} | {r['q_relative_to_global_min']:.3f} | {r['expected_signal_stars']:.2f} | {r['rate_at_grid_boundary']} |")
    lines+=['','The range of scenario optima is a sensitivity summary, not a confidence interval. Grid boundaries can hide larger shifts.',
        'Each curve below is re-zeroed at its own scenario minimum; use the table for between-scenario likelihood differences.','',
        '![Population profiles](population_profiles.png)','','![Visibility conversion](population_visibility.png)','',
        'Selected stars per escaped mass includes the orbit, field, time window and selection. It is not solely a luminosity-function correction.',
        'Candidate residual maps are descriptive calibration-mean diagnostics. They have no calibrated significance and do not replace the mixture likelihood.',
        'No q reweighting, literature-rate forcing, external population constraint or new confidence interval has been applied.',
        'Next inputs: escaping-population luminosity/mass-function constraints, foreground/control-region validation and calibrated recovery.',
        'Core mass functions require an explicit transfer model before being used as constraints on escaping stars.','']
    for c in comparisons:lines+= [f"![Candidate residuals]({c['artifact']}_residuals.png)",'']
    (out/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')
    assessment_path=grid/'production_assessment.json'
    if assessment_path.exists():
        previous=json.loads(assessment_path.read_text(encoding='utf-8'))
        (out/'previous_production_assessment.json').write_text(json.dumps(previous,indent=2),encoding='utf-8')
        previous.update(primary_rate_msun_myr=None,physical_rate_measurement_validated=False,
            headline='Absolute escape rate not externally calibrated; report conditional population scenarios.',
            population_sensitivity=population_summary(score['candidates']),population_diagnostics_report=str(out/'REPORT.md'))
        temp=assessment_path.with_name(assessment_path.name+'.population_review.tmp')
        temp.write_text(json.dumps(previous,indent=2),encoding='utf-8');temp.replace(assessment_path)
    (out/'population_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print('Population diagnostics complete:',out,flush=True);return summary


def watch(root,hours=72,interval=30):
    root=Path(root).resolve();deadline=time.monotonic()+hours*3600;attempts={}
    while time.monotonic()<deadline:
        for score in root.glob('*/production_grid_shared_v1/scores/joint_grid.json'):
            grid=score.parent.parent
            if not (grid/'production_assessment.json').exists() or (grid/'population_diagnostics_v1').exists():continue
            key=str(grid)
            if key in attempts:continue
            try:run(grid,selection=grid.parent/'selection')
            except FileExistsError:continue
            except Exception as exc:
                attempts[key]=str(exc);print('Population diagnostics failed:',grid,str(exc),flush=True)
                (grid/'population_diagnostics_error.json').write_text(json.dumps(dict(error=str(exc),source_sha256=sha256_file(Path(__file__))),indent=2),encoding='utf-8')
        time.sleep(min(interval,max(0,deadline-time.monotonic())))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--grid',type=Path);p.add_argument('--out',type=Path);p.add_argument('--selection',type=Path)
    p.add_argument('--survey-root',type=Path);p.add_argument('--hours',type=float,default=72);p.add_argument('--interval',type=float,default=30)
    a=p.parse_args()
    if a.survey_root:
        if min(a.hours,a.interval)<=0:p.error('Positive watch duration and interval required')
        watch(a.survey_root,a.hours,a.interval)
    elif a.grid:run(a.grid,a.out,a.selection)
    else:p.error('Specify --grid or --survey-root')
