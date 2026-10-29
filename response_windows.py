"""Build short lookback templates from the same full-response tracer contributions."""
import json
from pathlib import Path
import numpy as np
from provenance import sha256_file


def rebin(contributions,pixels,ages,window,n_pixels):
    c=np.asarray(contributions,float); p=np.asarray(pixels,int); ages=np.asarray(ages,float)
    if c.ndim!=2 or c.shape[1]!=len(p) or ages.shape!=p.shape or not np.isfinite(c).all() or (c<0).any():
        raise ValueError('Invalid tracer response')
    keep=(p>=0)&(p<n_pixels)&(ages<=window)
    return np.array([np.bincount(p[keep],weights=v[keep],minlength=n_pixels) for v in c])


def derive(full,out,duration):
    full,out=Path(full),Path(out)
    if out.exists(): raise ValueError('Use a fresh response-window directory')
    meta=json.loads((full/'response.json').read_text())
    maximum=meta['lookback_window']['duration_myr']
    if not 0<duration<maximum: raise ValueError('Need a shorter positive window')
    with np.load(full/'calibrated_counts.npz',allow_pickle=False) as data:
        counts={k:data[k] for k in data.files}
    with np.load(full/'tracer_response.npz',allow_pickle=False) as data:
        tracers={k:data[k] for k in data.files}
    desc=json.loads(str(counts['selection_description'].item()))
    history=desc.get('dynamics_mass_history') or {}
    if history.get('mode')!='prescribed_linear_escape': raise ValueError('Only prescribed linear histories supported')
    runpath=next((Path(k) for k in meta['input_sha256'] if Path(k).name=='run.json'),None)
    if runpath is None: raise ValueError('No dynamics provenance')
    run=json.loads(runpath.read_text())
    import pandas as pd
    from joint_grid import window_weights
    weights,window=window_weights(pd.DataFrame(dict(lookback_release_myr=tracers['lookback_myr'])),run,duration,tracers['mass_weight_msun'])
    fullstates=rebin(tracers['calibration_counts'],tracers['fit_pixel_index'],tracers['lookback_myr'],maximum,len(counts['counts']))
    if not np.allclose(fullstates,counts['signal_draws'],rtol=1e-10,atol=1e-10):
        raise ValueError('Full tracer contributions do not reconstruct full template')
    states=rebin(tracers['calibration_counts'],tracers['fit_pixel_index'],tracers['lookback_myr'],duration,len(counts['counts']))
    desc['lookback_window']=window
    counts.update(signal_draws=states,signal=states.mean(axis=0),selection_description=json.dumps(desc))
    with np.load(full/'response_maps.npz',allow_pickle=False) as d: maps={k:d[k] for k in d.files}
    maps['signal']=np.zeros_like(maps['signal']); maps['signal'][maps['fitmask']]=counts['signal']
    maps['signal_draws']=np.zeros((len(states),*maps['signal'].shape)); maps['signal_draws'][:,maps['fitmask']]=states
    # Full per-observation contributions were not stored. Do not claim a window MC error estimate.
    maps['mc_standard_error']=np.full_like(maps['signal'],np.nan)
    old=tracers['lookback_myr']>duration
    tracers['calibration_counts'][:,old]=0; tracers['expected_selected_stars'][old]=0
    tracers['mass_weight_msun']=weights
    out.mkdir(parents=True)
    np.savez_compressed(out/'calibrated_counts.npz',**counts)
    np.savez_compressed(out/'response_maps.npz',**maps)
    np.savez_compressed(out/'tracer_response.npz',**tracers)
    edges=np.linspace(0,run['config']['duration_myr'],21)
    ages=np.histogram(tracers['lookback_myr'],bins=edges,weights=tracers['expected_selected_stars'])[0]
    pd.DataFrame(dict(lookback_start_myr=edges[:-1],lookback_end_myr=edges[1:],expected_selected_stars=ages,
        cumulative_older_stars=np.cumsum(ages[::-1])[::-1])).to_csv(out/'age_contributions.csv',index=False)
    meta.update(lookback_window=window,expected_stars_in_fit_region=float(states.mean(axis=0).sum()),
        window_derivation='Exact sum of stored full-response fitted-region tracer contributions; no observation re-draw',
        cutflow_expected_counts=None,cache_hit=False,cache_key=None,
        observation_mc_error_available=False)
    for p in (full/'response.json',full/'tracer_response.npz'):
        meta['input_sha256'][str(p.resolve())]=sha256_file(p)
    meta['source_sha256']['response_windows.py']=sha256_file(Path(__file__))
    meta.update(template_sha256=sha256_file(out/'calibrated_counts.npz'),maps_sha256=sha256_file(out/'response_maps.npz'))
    (out/'response.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    return meta
