# 安裝與資料準備

本儲存庫只包含研究程式、測試、設定、少量合成範例與驗證摘要。公開程式不等於已驗證的物理 MLR 測量。

## 安裝

已於 Python 3.11 測試。

```powershell
git clone https://github.com/yohan4567/Tidal-Tail-Analysis.git
cd Tidal-Tail-Analysis
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

Linux/macOS 改用 `.venv/bin/python`。galpy 原生 C extension 無法載入時，本專案有較慢的 CPU 積分路徑。

## 不需要外部資料的核心測試

```powershell
New-Item -ItemType Directory -Force outputs | Out-Null
.\.venv\Scripts\python.exe -m unittest test_model test_extent test_dynamics test_mass_evolution test_joint_grid test_profile_recovery test_dense_convergence test_validation
```

這些是軟體與數值測試，不需下載 Gaia 場或執行完整星流實驗。Linux/macOS 先執行 `mkdir -p outputs`。

## 外部選星程式

觀測串接重用原研究的 `tidal_study_v2`，沒有把另一個專案的原碼、私人 repo 或觀測資料複製到此儲存庫。請自行取得有權使用的版本，設定其路徑：

```powershell
$env:GC_TIDAL_STUDY = 'D:\research\tidal_study_v2'
```

Linux/macOS：

```bash
export GC_TIDAL_STUDY=/path/to/tidal_study_v2
```

該目錄至少需包含 `step3_tail_search.py` 與原程式匯入的相依檔案。實際場星觀測流程還會使用其 `catalogs/`、`gaia_fields/` 及對應來源中繼資料。支援的個別 CLI 也可傳 `--study`。舊本機資料夾排列仍保留為 fallback，但不是新安裝的必要排列。

外部程式與資料可用後，再跑完整整合測試：

```powershell
.\.venv\Scripts\python.exe -m unittest test_profile_recovery test_dense_convergence test_joint_grid test_mass_evolution test_validation test_science_calibration test_observation test_model test_selection test_extent test_dynamics
```

文獻與本機 60 項測試的結果見 [DENSE_TEST_RESULTS.md](DENSE_TEST_RESULTS.md)。該數字包含外部選星整合測試，不能宣稱乾淨 clone 不需資料就能跑全部整合功能。

## 資料與校正

1. 從原始來源準備星團目錄、Gaia 場及中繼資料；來源與查詢流程見 [SCIENCE_CALIBRATION.md](SCIENCE_CALIBRATION.md)。
2. 下載官方 Gaia 完整率地圖，見 [OBSERVATION_V2.md](OBSERVATION_V2.md)。
3. 準備有來源的等時線、人口質量定義、品質校正與誤差母樣本。不要用 `examples/*_SYNTHETIC.csv` 報告真實 MLR。
4. 將可重新產生的資料與結果放在 `outputs/`。此目錄被 Git 忽略。

`configs/Pal5.json` 是採用的星團參數快照，不是觀測擬合結果。完整 joint-grid、驗證與後續測試命令見 [JOINT_GRID.md](JOINT_GRID.md)；其中指向過去 `outputs/...` 的命令需要先生成那些輸入。

## 公開內容與排除內容

保留：Python 程式、測試、依賴清單、靜態設定、方法文件、合成小範例，以及 `docs/results/`、`docs/figures/` 中經挑選的彙總結果。

排除：`.env`、憑證、虛擬環境、Python 快取、原始 Gaia CSV/FITS、完整粒子與 NPZ、下載快取、暫存檔、log、完整 outputs 及其他專案副本。沒有刪除本機研究資料。

本機完整來源雜湊與軌跡仍留在 outputs；公開的圖與小表是結果摘要，不是所有實驗輸入的替代品。
