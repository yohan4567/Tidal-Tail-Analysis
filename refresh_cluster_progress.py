"""Refresh file-based survey progress without modifying or restarting scientific jobs."""
import argparse,json,re,time
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
from run_cluster_survey import SAMPLE


def refresh(root):
    root=Path(root).resolve();rows=[]
    for display,cid in SAMPLE:
        dest=root/cid;audit=dest/'input_audit.json'
        entry=json.loads(audit.read_text()) if audit.exists() else dict(cluster=cid,display_name=display,status='queued')
        link=dest/'reuse.json'
        if link.exists():entry['reused_outputs']=json.loads(link.read_text())
        rows.append(entry)
    report=dict(clusters=len(SAMPLE),entries=rows,updated_utc=datetime.now(timezone.utc).isoformat())
    for entry in report['entries']:
        dest=root/entry['cluster'];links=entry.get('reused_outputs',{})
        source=links.get('production_grid_source')
        if source:
            grid=Path(source);plan_path=grid/'plan.json'
            if plan_path.exists():
                plan=json.loads(plan_path.read_text());jobs=[];integrated=0;complete=0
                particles=2*round(max(plan['times_myr'])/plan['release_interval_myr'])
                for rate in plan['rates']:
                    for seed in plan['seeds']:
                        run=grid/f'rate_{rate:g}_seed_{seed}';done=(run/'simulation/run.json').exists()
                        count=particles if done else 0;log=run/'simulation.log';stamp=None
                        if log.exists():
                            data=log.read_text(encoding='utf-8',errors='replace')
                            matches=re.findall(r'Integrated (\d+)/(\d+) tracers',data)
                            if matches and not done:count=min(particles,int(matches[-1][0]))
                            stamp=datetime.fromtimestamp(log.stat().st_mtime,timezone.utc).isoformat()
                        if done:complete+=1
                        integrated+=count
                        if log.exists() or done:jobs.append(dict(rate=rate,seed=seed,complete=done,last_reported_tracers=count,log_updated_utc=stamp))
                total=len(plan['rates'])*len(plan['seeds'])
                entry['dynamics_progress']=dict(completed_runs=complete,total_runs=total,last_reported_tracers=integrated,
                    total_tracers=total*particles,tracer_fraction=integrated/(total*particles),jobs=jobs,
                    grid_scores_available=(grid/'scores/joint_grid.json').exists(),
                    caveat='Tracer count fraction is not elapsed-time or overall scientific-validation completion. Log output may be buffered.')
        old=root.parent/'nine_cluster_real_inputs_v1'/entry['cluster']
        state=old/'parent_download_status.json'
        if state.exists():entry['parent_download']=json.loads(state.read_text())
        parent=links.get('calibration_parent_source')
        if parent:entry['calibration_parent_source']=parent
    report['status_counts']=dict(Counter(e['status'] for e in report['entries']))
    report['refresh_utc']=datetime.now(timezone.utc).isoformat()
    def write(path,text):
        temp=path.with_name(path.name+'.refresh.tmp');temp.write_text(text,encoding='utf-8');temp.replace(path)
    write(root/'progress.json',json.dumps(report,indent=2))
    labels={'input_preparation_failed':'缺少輸入或準備失敗','observations_ready_calibration_pending':'選星完成，等待自身校正','conditional_grid_inputs_available':'條件式網格輸入齊備'}
    lines=['# 29 個星團執行進度','',f'最新刷新 UTC：{report["refresh_utc"]}','',
        '此表讀取實際輸出及日誌；準備完成不等於完成 MLR 擬合。','',
        '| 星團 | 輸入與選星階段 | 擬合區星數 | 動力學與校正進度 |','|---|---|---|---|']
    for e in report['entries']:
        detail=[]
        if e.get('dynamics_progress'):
            d=e['dynamics_progress'];detail.append(f"已完成 {d['completed_runs']}/{d['total_runs']} 組；已回報 {d['last_reported_tracers']:,}/{d['total_tracers']:,} 粒子（{100*d['tracer_fraction']:.2f}%）")
        if e.get('parent_download'):detail.append('校正母樣本下載：'+e['parent_download']['status'])
        if e.get('error'):detail.append(e['error'].split(':')[0])
        lines.append(f"| {e['display_name']} | {labels.get(e['status'],e['status'])} | {e.get('n_fit_stars','—')} | {'；'.join(detail) or '尚未開始 MLR 擬合'} |")
    lines+=['','粒子比例不是總耗時或科學驗證完成度。緩衝輸出可能使數字暫時停留；程序是否存活需另行確認。',
        'Pal 5 沿用原執行目錄，來源見 Pal5/reuse.json。缺少人口、觀測校正或完整率驗證的星團不會自動啟動正式 MLR 推論。']
    write(root/'PROGRESS.md','\n'.join(lines)+'\n')
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--watch',action='store_true');p.add_argument('--interval',type=float,default=30)
    p.add_argument('--hours',type=float,default=72)
    a=p.parse_args()
    if a.interval<5 or a.hours<=0:p.error('Use interval >= 5 seconds and positive hours')
    end=time.monotonic()+a.hours*3600
    while True:
        result=refresh(a.root);print(result['refresh_utc'],result['status_counts'],flush=True)
        if not a.watch or time.monotonic()>=end:break
        time.sleep(a.interval)

if __name__=='__main__':main()
