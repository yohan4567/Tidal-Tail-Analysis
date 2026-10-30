"""Prepare and report the 29-cluster sample without duplicating live simulations."""
import json,shutil,sys
from pathlib import Path
from datetime import datetime,timezone
import prepare_cluster_batch as prep
HERE=Path(__file__).resolve().parent
SAMPLE=[('ESO 452-SC11','ESO452-SC11'),('NGC 288','NGC288'),('NGC 362','NGC362'),('NGC 1851','NGC1851'),
('NGC 2298','NGC2298'),('NGC 3201','NGC3201'),('NGC 5466','NGC5466'),('M80','NGC6093'),('M4','NGC6121'),
('NGC 6139','NGC6139'),('M92','NGC6341'),('NGC 6366','NGC6366'),('NGC 6397','NGC6397'),('M30','NGC7099'),
('Pal 5','Pal5'),('M5','NGC5904'),('NGC 6144','NGC6144'),('M107','NGC6171'),('NGC 6535','NGC6535'),
('Omega Cen','NGC5139'),('NGC 5286','NGC5286'),('M62','NGC6266'),('M19','NGC6273'),('NGC 6304','NGC6304'),
('NGC 6355','NGC6355'),('NGC 6541','NGC6541'),('M22','NGC6656'),('M15','NGC7078'),('M2','NGC7089')]

def progress(root):
    rows=[]
    for display,cid in SAMPLE:
        dest=root/cid;audit=dest/'input_audit.json'
        row=json.loads(audit.read_text()) if audit.exists() else dict(cluster=cid,display_name=display,status='queued')
        link=dest/'reuse.json'
        if link.exists():row['reused_outputs']=json.loads(link.read_text())
        rows.append(row)
    report=dict(updated_utc=datetime.now(timezone.utc).isoformat(),clusters=len(SAMPLE),entries=rows,
        note='Prepared real observations are not completed MLR inference; unavailable calibrations are explicit.')
    (root/'progress.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    lines=['# 29 個星團執行進度','',f'更新時間 UTC：{report["updated_utc"]}','',
           '| 星團 | 識別碼 | 輸入與選星狀態 | 擬合區星數 |','|---|---|---|---|']
    for r in rows:lines.append(f'| {r["display_name"]} | {r["cluster"]} | {r["status"]} | {r.get("n_fit_stars","—")} |')
    lines+=['','Pal 5 的計算來源見其 reuse.json；不重複積分。其他星團只有在自身人口、誤差及品質校正齊備後才進行 MLR 擬合。',
            '原始形態標籤若另行匯入，僅作啟發式候選分類，不等於校準後的尾巴偵測顯著性。']
    (root/'PROGRESS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return report


def main():
    root=HERE/'outputs/cluster_survey_29_v1';root.mkdir(parents=True,exist_ok=True)
    old=HERE/'outputs/nine_cluster_real_inputs_v1'
    for display,cid in SAMPLE:
        dest=root/cid;dest.mkdir(exist_ok=True);src=old/cid
        reused={}
        if (src/'selection/selection.json').exists() and not (dest/'selection').exists():
            shutil.copytree(src/'selection',dest/'selection');reused['selection_source']=str(src/'selection')
        if (src/'adopted_dynamics.json').exists() and not (dest/'adopted_dynamics.json').exists():
            shutil.copy2(src/'adopted_dynamics.json',dest/'adopted_dynamics.json')
        if (src/'parent_random10.csv.gz').exists():reused['calibration_parent_source']=str(src/'parent_random10.csv.gz')
        if cid=='Pal5':reused['production_grid_source']=str(src/'production_grid')
        if reused:(dest/'reuse.json').write_text(json.dumps(reused,indent=2),encoding='utf-8')
    progress(root)
    for item in SAMPLE:
        prep.CLUSTERS=[item]
        prep.prepare(root,workers=1)
        progress(root)
    import pandas as pd
    master=pd.read_csv(prep.DEFAULT_STUDY/'catalogs/master_clusters.csv').set_index('id')
    master.loc[[cid for _,cid in SAMPLE],['ra','dec','R_sun','pmra','pmdec','rv','mass','r_hm_pc','feh','ebv']].to_csv(root/'adopted_parameters.csv')
    full=progress(root)
    (root/'batch_audit.json').write_text(json.dumps(full,indent=2),encoding='utf-8')
    print('Batch observation preparation complete; see PROGRESS.md',flush=True)

if __name__=='__main__':main()
