"""Read-only audit of inference outputs: computation is distinct from mass closure."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path

import numpy as np
from provenance import sha256_file
from population_model import Population


def audit_posterior(path):
    path=Path(path).resolve();report=json.loads(path.read_text(encoding='utf-8'))
    prior=report.get('population_prior') or {};constraint=prior.get('constraint',{})
    measured_tail=constraint.get('measurement_scope')=='present_day_observed_tail_segment'
    support=prior.get('support_audits',[])
    population_checks=[]
    for config in prior.get('response_configs',[]):
        cfg=Path(config)
        try:
            spec=json.loads(cfg.read_text(encoding='utf-8'))
            folder=(cfg.parent/spec['population']).resolve();pop=Population.load(folder)
            total=float(np.sum(pop.table.current_mass_msun*pop.table.number_per_msun))
            population_checks.append(dict(config=str(cfg),metadata_sha256=sha256_file(folder/'population.json'),
                table_hash_verified=True,current_mass_normalization=total,
                unit_mass_pass=bool(np.isclose(total,1.,atol=1e-10,rtol=0)),
                synthetic=bool(pop.metadata.get('synthetic_calibration',False))))
            if not prior.get('support_audits') and constraint.get('measured_mass_range_msun'):
                measured=constraint['measured_mass_range_msun'];represented=pop.metadata.get('represented_ms_mass_range_msun')
                if represented:
                    support.append(dict(mf_extrapolated=represented[0]<measured[0] or represented[1]>measured[1],
                        measured_mass_range_msun=measured,represented_ms_mass_range_msun=represented))
        except (OSError,ValueError,KeyError) as exc:
            population_checks.append(dict(config=str(cfg),error=str(exc),unit_mass_pass=False))
    ranges=[]
    for r in report['results'].values():
        finite=all(np.isfinite([r[a]['median'],*r[a]['credible68'],*r[a]['credible95']]).all()
                   for a in ('rate','lookback_time'))
        ranges.append(dict(finite=bool(finite),rate=r['support']['rate_msun_myr'],time=r['support']['time_myr'],
            rate_edge_quadrature_mass=[r['rate']['first_node_quadrature_mass'],r['rate']['last_node_quadrature_mass']],
            time_edge_quadrature_mass=[r['lookback_time']['first_node_quadrature_mass'],r['lookback_time']['last_node_quadrature_mass']]))
    # Existing metadata explicitly declares the unresolved scope; never infer that
    # f=1, non-synthetic inputs or passing software tests establish physical closure.
    checks=[
        dict(stage='population_unit_mass',status=('pass' if all(p['unit_mass_pass'] for p in population_checks) else 'fail') if population_checks else 'not_checked',
             evidence='Population CSV checksum and numerical unit living-mass normalization; not total physical mass closure'),
        dict(stage='posterior_numerics',status='pass' if all(r['finite'] for r in ranges) else 'fail',
             evidence='Finite medians and intervals only; not a recovery test'),
        dict(stage='escaping_population',status='partial_constraint' if measured_tail else 'unresolved',
             evidence='Measured tail segment; whole-field/time transfer still assumed' if measured_tail else 'Global MF plus an uncalibrated escape-transfer scenario'),
        dict(stage='unobserved_mass',status='unresolved',evidence='No externally supported low-mass/remnant mass closure in this pipeline'),
        dict(stage='dynamical_mass_conversion',status='unresolved',evidence='Living-star/dynamical conversion remains prescribed, not inferred'),
        dict(stage='selection_and_background',status='validation_required',evidence='Per-star/field calibration does not establish transfer to the escaping population or background shape'),
        dict(stage='zero_signal_model',status=(report.get('background_only_model') or {}).get('status','not_computed'),
             evidence='A positive-only MLR prior cannot establish tail detection'),
        dict(stage='independent_recovery_and_convergence',status='validation_required',evidence='Software tests and fixed-template recovery are insufficient for model-mismatch validation'),
        dict(stage='dynamics_and_footprint',status='validation_required',evidence='Orbit/halo/scale/history applicability and finite survey footprint remain conditional')]
    return dict(posterior=str(path),posterior_sha256=sha256_file(path),checks=checks,
        posterior_support=ranges,measured_tail_constraint=measured_tail,
        population_input_checks=population_checks,
        mf_extrapolation_detected=any(a.get('mf_extrapolated',False) for a in support) if support else None,
        primary_total_absolute_rate_msun_myr=None,total_absolute_rate_supported=False,
        allowed_output='Conditional rate distribution for the declared represented living-star model',
        interpretation='Audit of available evidence, not an independent physical validation or a proof of zero mass loss.')


def run(root,out):
    root=Path(root).resolve();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    reports=[audit_posterior(p) for p in sorted(root.glob('*/empirical_mf_v1/*/posterior/posterior_summary.json'))]
    reports += [audit_posterior(p) for p in sorted(root.glob('*/tail_mf_literature_v1/*/posterior/posterior_summary.json'))]
    reports += [audit_posterior(p) for p in sorted(root.glob('*/empirical_direct_v1/*/posterior/posterior_summary.json'))]
    result=dict(updated_utc=datetime.now(timezone.utc).isoformat(),reports=reports,
        source_sha256=sha256_file(Path(__file__)),scope='Existing completed results only; no simulation restarted or relabelled')
    (out/'audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    lines=['# 實驗流程稽核','','完成計算不等於完成絕對質量閉合。下表僅檢查既有結果的證據與限制。','',
           '| 結果 | 逃逸人口 | 外推被標示 | 絕對總流失率 |','|---|---|---|---|']
    for r in reports:
        name=str(Path(r['posterior']).relative_to(root))
        extrap='未能檢查' if r['mf_extrapolation_detected'] is None else ('有外推' if r['mf_extrapolation_detected'] else '未檢出外推')
        lines.append(f"| {name} | {'尾巴測光；轉移仍是假設' if r['measured_tail_constraint'] else '全域 MF＋δ 假設'} | {extrap} | 尚無充分支持 |")
    lines+=['','完整檢查、人口單位質量／CSV 檢查碼及摘要檔雜湊見 audit.json。缺少範圍資料會標為未能檢查，不視為通過。',
            '缺少資料不能由改變標籤、縮小誤差或增加運算次數補足。']
    (out/'AUDIT_ZH.md').write_text('\n'.join(lines),encoding='utf-8')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=run(a.root,a.out)
    print(json.dumps(dict(completed_posteriors=len(r['reports']),total_absolute_rate_supported=False)))
