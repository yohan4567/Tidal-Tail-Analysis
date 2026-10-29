# Tidal Tail Analysis

球狀星團潮汐尾的 particle-spray 模擬、Gaia 觀測模型與條件式質量流失率推論。

**新安裝請先閱讀 [安裝與資料準備](docs/SETUP_AND_DATA.md)。** 此儲存庫保留程式、測試、設定與驗證摘要；Gaia 資料、外部選星專案、快取及完整模擬輸出需另行準備。過去實驗命令中的 `outputs/...` 是本機研究產物，不是 clone 後已存在的檔案。

T–MLR 聯合流程的實測、通過與未通過項目，見 [實際驗證結果](docs/JOINT_VALIDATION_RESULTS.md)。目前總計數穩定，但細空間模板與區間校準尚未全面通過；最佳值仍屬條件式試驗。

文獻對照、分箱診斷與改善優先順序，見 [改善研究](docs/IMPROVEMENT_RESEARCH_ZH.md)。其中另澄清完整候選點涵蓋率與 MLR profile 區間的差別。

改善後的實測結果見 [8000 粒子與 Profile 校準測試](docs/DENSE_TEST_RESULTS.md)：已執行三種子密集模型、4 Gyr 稀疏年代預演及獨立擴大校準。仍需區分條件式數值通過與物理模型驗證。

球狀星團 particle-spray 模擬、原專案選星、尾巴尺度分析與條件式流失率擬合。
這是一個可執行的研究工具第一版；真實流失率仍需要外部人口與完整度校正。

**目前推薦的觀測路徑是 Hybrid v2**：已加入人口模型、官方 TCG 完整度地圖、缺測／品質機率、相關誤差、加權模板與多模型比較。
完整操作見 [Observation v2](docs/OBSERVATION_V2.md)。下方 `calibrated` 模式保留為舊版 donor 近似的對照；新路徑使用 `--mode weighted`。

正式人口輸入與實際 Gaia 母樣本校正、獨立保留樣本檢查及敏感度模型，見 [科學校正流程](docs/SCIENCE_CALIBRATION.md)。`examples/create_v2_demo.py` 仍僅供軟體測試；不會自動升格為正式校正。

九組人口／零點敏感度、粒子數與觀測抽樣比較、重新擬合的 Poisson bootstrap，見 [科學驗證](docs/VALIDATION.md)。這些診斷會區分數值穩定性、假設敏感度與條件式擬合檢查。

原始設計與尚需校正的研究項目見 [詳細製作規格](docs/IMPLEMENTATION_PLAN.md)；已實作狀態以 v2 操作文件為準。

## 快速開始

在 PowerShell 進入本資料夾：

```powershell
cd "path\to\Tidal-Tail-Analysis"
.\.venv\Scripts\python.exe run_pipeline.py --cluster Pal5 --legacy --out outputs/Pal5_analysis
```

這個預設模式不需要模擬或人口校正，會完成：

1. 直接重用原專案的 Gaia 選星與遮罩。
2. 輸出原始計數、背景圖與預先指定的雙向軌道。
3. 沿軌道測量條件式尾巴候選段、可追蹤延伸距離與橫向寬度。
4. 產生 `REPORT.md`、`manifest.json`、圖表與 CSV。

每次使用新的輸出資料夾；程式拒絕覆寫舊結果。失敗時保留 manifest 與各階段 log。
預設不會輸出實測 mass-loss rate。

## 安裝

本機驗證 Python 3.11。已有 `.venv` 時不用重建。

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

若 ensurepip 失敗，可執行 `python -m venv --without-pip .venv`，
再用有 pip 的 Python 執行 `python -m pip --python .venv/Scripts/python.exe install -r requirements.txt`。

galpy C extension 在本機缺少可載入的相依庫；程式自動使用 CPU 路徑，
釋放抽樣由 galpy 完成，粒子由快取軌跡的 Cartesian SciPy DOP853 積分，不使用 GPU。
C extension 警告本身不代表模擬失敗；原生 C 路徑保留，但未在本機驗證。
從本資料夾執行，讓 `.galpyrc` 關閉不必要的網路版本查詢。

## 依賴你原本的選星專案

預設讀取相鄰的：

`../TidaleffectsOnGCs-main/TidaleffectsOnGCs-main/tidal_study_v2/`

可用 `--study` 指定另一個路徑。模組直接呼叫 `step3_tail_search.py` 的函式，
不另外重寫 RUWE、flux excess、視差、自行、CMD、鄰近物體與缺失 tiles 遮罩。
CMD ridge 與軌道修正是否啟用均由真實資料決定，再固定套用於模擬觀測。
輸出會記錄來源資料、程式 SHA256 與每一步通過的星數。

`--legacy` 明確使用原始天區，不使用 v3 規畫。
本機 Pal5 v3 計畫要求的 10 度資料不存在，原始 5 度資料存在，所以範例使用此選項。
程式不會默默換成舊資料。舊 Pal5 檔有 39,677 顆星，原程式選後 2,117 顆，
新介面已逐列核對一致。其 PM 下載窗口較窄，原程式會使用固定自行切選，原因保留在報告。

## 模擬星團與模型尺度

```powershell
.\.venv\Scripts\python.exe run_pipeline.py --cluster Pal5 --legacy --simulate --duration-myr 100 --particles 20 --out outputs/Pal5_smoke
```

20 粒子只用於確認流程。正式研究需要更多粒子、不同 seed、不同模擬時間的收斂檢查。
也可單獨執行：

```powershell
.\.venv\Scripts\python.exe stream_model.py --config configs/Pal5.json --duration-myr 1000 --particles 2000 --out outputs/Pal5_1gyr
.\.venv\Scripts\python.exe model_extent.py --run outputs/Pal5_1gyr --out outputs/Pal5_model_extent
```

模型假設：

- 使用 galpy 1.12 的 `chen24spraydf`，兩臂各一半 tracer，均勻抽樣過去釋放時間。
- 銀河勢為 MWPotential2014，保留原生尺度 ro=8 kpc、vo=220 km/s。
- 星團為移動的固定質量 Plummer 勢；`a = r_hm * sqrt(2^(2/3)-1)`。
- galpy 觀測轉換採 z_sun=20.8 pc、特殊速度 [-11.1,12.24,7.25] km/s。
- 單粒子質量權重為 `mdot * duration / n_particles`；mdot 是兩臂總和，沒有額外乘二。
- 改 mdot 只改權重；目前不自洽演化星團質量、尺度或逃逸歷史，也沒有粒子間引力與內部鬆弛。
- 若累積流失質量不可忽略，需要演化星團勢。run.json 記錄代表質量占固定星團質量的比例。
- 原選星程式的 skymath 保留其自己的座標框架與軌道設定；模擬與選星的框架假設須納入敏感度研究。

`model_extent.json` 回報每臂 95% 質量的角向徑向尺度與模型 3D 距離。
這是從中心量的分位數，不是沿尾巴長度、觀測偵測距離或 Jacobi 半徑。
少於 100 粒子／臂會標示 small_sample；超過也不代表已完成收斂。

模擬輸出：`particles.csv`、`orbit.csv`、`overview.png`、`run.json`。
CSV 含今天的位置、距離、自行、RV、釋放回溯時間、質量權重與模型臂標記。
軌道時間負值代表過去，Cartesian 軸使用 galpy 慣例。

## 觀測尾巴的長度與寬度

```powershell
.\.venv\Scripts\python.exe tail_extent.py --selection outputs/Pal5_analysis/selection --out outputs/Pal5_profile --bin-deg 0.4 --half-width-deg 0.3
```

方法與回報方式：

- 軌道事先由 catalogue 六維資訊產生，不從同一張密度圖挑最亮的方向。
- 在 gnomonic 投影沿軌道切互不重疊的 bins，計算原始星數與原程式背景。
- 預設排除中心 1.5*r_t，搜尋至 0.9*min(field radius,5 deg)。
- 預設 bin 寬 max(0.25 deg,3 pixels)，橫向半寬 max(0.25 deg,2.1 pixels)，均可改。
- 每個 bin 以固定已知背景的 Poisson 上尾機率比較，對兩方向全部有效 bins 作 Bonferroni 校正。
- 預設至少 70% 條帶面積有效，且須連續兩個候選 bins。遮罩或軌道／視野終止會打斷連通距離。
- 回報從內邊界連續可追蹤的 extent，以及另外的分離候選段。缺乏訊號時輸出 null，不給 0 長度。
- extent bracket 只是端點所在 bin 的解析度範圍，不是信賴區間。
- 寬度是候選 bin 內扣背景後的橫向第二矩 RMS，受條帶截斷與像素影響；不是去卷積 Gaussian sigma。

forward/backward 是軌道方向，不直接宣稱為真實 leading/trailing。
角尺度是投影平面弧長；pc 換算採星團距離的 D*theta，不是 3D 去投影。
無連通候選段不代表沒有潮汐尾。此方法也不是通用的 halo、shell 或非軌道對齊結構搜尋器。

**統計限制：** 背景由同一資料估計，可能吸收尾巴；目前未傳遞背景形狀誤差。
因此 p 值與候選段只是在固定背景與軌道假設下的條件式結果。
校正不涵蓋掃描多個方向、寬度、設定或星團的多重比較，軌道式選星本身也可能造成偏好。
圖上的 sqrt(N+B) 誤差棒只是探索用途，不能視為完整背景不確定度。

## 注入恢復與假訊號檢查

```powershell
.\.venv\Scripts\python.exe validate_tail.py --selection outputs/Pal5_analysis/selection --true-extent-deg 2 --selected-stars 300 --sigma-deg 0.1 --trials 100 --out outputs/Pal5_recovery
```

此處 300 是選後、未遮罩天區內的兩臂總期望星數，不是質量或流失率。
實驗在固定背景圖上產生 Poisson 純背景與注入訊號，輸出候選比例、二項比例區間、
恢復端點與偏差。這只是 selected-space 測試，不會重建 CMD、重擬背景或驗證 Gaia 完整度。
請針對強度、寬度、長度及設定建立網格，而非用單一強注入宣稱整體完整度。

本機一次固定 seed 的測試：100 次純背景未出現候選段，95% 比例區間上端約 3.6%；
注入 300 顆、2 度長度的測試 100 次均找到連通段，端點偏差約在一個 bin 內。
這不是實際天空的假訊號率保證。

## 有外部校正時的流失率流程

需要兩份 CSV：

| 檔案 | 欄位與意義 |
|---|---|
| population.csv | M_G、本徵 bp_rp、number_per_msun：每個光度格的恆星數／總逃逸質量 |
| completeness.csv | apparent G、completeness：科學品質切選之前的偵測／入表機率 |

人口模型應對應逃逸族群；分母須考慮看不到的星與殘骸。
completeness 必須涵蓋所有模擬星等，避免重複扣除後續 RUWE／CMD 等切選。
需要額外提供固定 A_G、E(BP-RP) 及來源說明；程式沒有內建真實 Pal5 的校正值。

```powershell
.\.venv\Scripts\python.exe run_pipeline.py --cluster Pal5 --legacy --mode calibrated --reuse-run outputs/Pal5_smoke/simulation --population population.csv --completeness completeness.csv --a-g 0.1 --reddening 0.05 --calibration-description "Describe population and completeness sources" --null-trials 200 --out outputs/Pal5_calibrated
```

消光數字只是命令語法範例，須換成實際校正。正式運算也須替換小樣本模擬。
`examples/*_SYNTHETIC.csv` 是任意軟體測試資料，使用時必須加 `--synthetic-calibration`，
它們不適用於回報物理流失率，結果報告會明確標記。

模擬觀測的條件式近似：

- 每 tracer 依恆星數抽一個光度格，轉成 apparent G 與顏色。
- expected_count = tracer mass * sum(number_per_msun) * completeness。
- 在原始 Gaia field 的鄰近 G／顏色中抽取誤差、RUWE、flux-excess 殘差組合。
- 對 PM 與視差加入獨立 Gaussian 誤差，再套用原專案的下載窗口、品質、CMD 與遮罩。
- 尚無誤差協方差、空間系統誤差、光度散布或差異消光；donor 也受原下載窗口限制。
- 增加 tracer 數並重複 seed，評估人口／誤差抽樣與模板 Monte Carlo 噪聲。

選後模板以預期星數加權；mock 列數不是實際星數，原始 tracer 質量也不是星數模板。

擬合模型為 `N_j ~ Poisson(A*S_j + beta*B_j)`，同時估 A 與非負背景倍率。
輸出 rate=A*reference_rate、68% profile likelihood 區間、單側 95% 上限。
區間使用漸近 likelihood 門檻，在零邊界與低計數時只為近似；不包括人口、完整度或模板形狀誤差。
`--null-trials` 可對固定模板與背景形狀作 parametric null bootstrap；不代表盲搜尋顯著性。
沒有選後訊號或模板不可辨識時，程式報錯，不製造流失率。

獨立合成計數測試：

```powershell
.\.venv\Scripts\python.exe fit_rate.py --demo --null-trials 100 --out outputs/synthetic_rate_demo.json
```

## 程式入口與驗證

| 程式 | 用途 |
|---|---|
| run_pipeline.py | 統一入口、階段 log、manifest、結果報告 |
| stream_model.py | Chen 模擬與數值軌道檢查 |
| make_config.py | 從原 catalogue 匯入另一顆星團 |
| observe_particles.py | 外部校正的模擬觀測 |
| selection_bridge.py | 重用原選星，輸出未平滑計數與模板 |
| tail_extent.py | 沿軌道的條件式長度／寬度 |
| model_extent.py | 模擬粒子的質量分位尺度 |
| validate_tail.py | 固定背景的選後注入恢復 |
| fit_rate.py | 條件式流失率、profile 區間與 null bootstrap |

```powershell
.\.venv\Scripts\python.exe -m unittest -v test_model test_selection test_extent
```

測試涵蓋兩臂質量正規化、觀測權重與零誤差極限、凍結 CMD／遮罩／下載窗口、
曲線投影、已知尾巴端點恢復、純背景、遮罩中斷、加權分位數與 rate 區間／零邊界。
實際資料選星另由 selection_bridge 每次逐列核對原程式。

正式科學結果還需：真實人口與完整度校正、背景錯配與空間誤差模型、
多組參數／seed 的收斂研究、完整觀測層的注入恢復，以及對銀河勢與星團演化的敏感度分析。

## 參考

- Chen et al. particle spray：https://arxiv.org/abs/2408.01496
- 作者範例：https://github.com/ybillchen/particle_spray
- 星流與質量流失研究：https://arxiv.org/abs/2411.19899

程式直接呼叫套件，不複製作者的實作。
# Joint time and mass loss experiments

The constant-rate, mass-history-coupled T–MLR workflow is documented in
[docs/JOINT_GRID.md](docs/JOINT_GRID.md), with the scientific protocol in
[docs/EXPERIMENTAL_METHOD_ZH.md](docs/EXPERIMENTAL_METHOD_ZH.md).
It supports stratified release sampling, Tmax window reuse, conditional profiles,
independent-template noise diagnostics and frozen-template injection recovery.
These tools do not automatically validate a physical mass-loss measurement.
