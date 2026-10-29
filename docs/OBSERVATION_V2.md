# Hybrid observation pipeline v2

本版完成人口模型、Gaia 空間完整度、缺測／聯合品質機率、相關天體測量誤差、選後期望計數、共享校正抽樣的擬合與模型比較接口。
程式可執行，不代表真實星團的人口／品質校正已經完成。只有合成校正的結果會持續標示 SYNTHETIC。

正式 PARSEC 下載、Gaia 母樣本與保留樣本驗證，以及逃逸質量分布／零點敏感度設定，見 [科學校正流程](SCIENCE_CALIBRATION.md)。人口匯入現在排除 PARSEC 的暗終止節點；這不會自動補出完整殘骸族群。

## 1. 安裝與 Gaia 地圖

目前本機 Windows 上 GaiaUnlimited 0.3.3 的 healpy 相依套件無法編譯。
因此推薦使用 `tcg_map`：直接讀官方 TCG M10 地圖並評估同一發表模型，使用 astropy-healpix，沒有偽裝成已安裝原生 GaiaUnlimited。
也保留 `gaiaunlimited` 後端供能安裝該套件的環境使用。

```powershell
python -m pip --python .venv/Scripts/python.exe install -r requirements.txt
.\.venv\Scripts\python.exe selection_model.py fetch-gaia-map --mode multi --out outputs/cache/gaia/allsky_uniq_10.fits
```

原生套件為可選安裝：`requirements-observation.txt`。
地圖約 36 MB，已下載時不重跑 fetch；程式拒絕覆寫。地圖與其 SHA256 sidecar 都在忽略版控的快取。
模型依據：[TCG 論文](https://arxiv.org/abs/2208.09335)、[官方介面](https://gaiaunlimited.readthedocs.io/en/latest/api/gaiaunlimited.selectionfunctions.DR3SelectionFunctionTCG.html)。

## 2. 可重現的軟體示範

以下人口、品質與誤差資料是任意合成設定；官方完整度地圖是真實參考資料。

```powershell
.\.venv\Scripts\python.exe examples/create_v2_demo.py --out outputs/demo_v2 --real-gaia-map outputs/cache/gaia/allsky_uniq_10.fits
.\.venv\Scripts\python.exe run_pipeline.py --cluster Pal5 --legacy --mode weighted --simulate --duration-myr 1000 --particles 200 --response-config outputs/demo_v2/response.json --out outputs/Pal5_demo_v2
```

可用 `--reuse-run outputs/Pal5_demo_v2/simulation` 重用同一動力學。
這是軟體串接示範，不能從該結果報告 Pal5 的物理流失率。原始 field 的 2,117 顆候選星仍逐列核對原程式。

## 3. 匯入真正的人口模型

`population_model.py` 支援 PARSEC 註解表頭的文字檔及 canonical CSV。
PARSEC 檔案必須選擇 Gaia 相容光度系統，而且只含一個年齡／化學組成；程式檢查 logAge/MH 與宣告參數是否一致。
輸入星等必須是未消光的絕對星等，避免在 response 階段再加一次消光。

```powershell
.\.venv\Scripts\python.exe population_model.py --isochrone isochrone.dat --age-years 12000000000 --metallicity-dex -1.4 --mass-basis isochrone_current_stellar_mass --source-description "Actual model version, abundance assumptions and source" --out outputs/population_calibrated
```

上述年齡與金屬量只示範語法，須換成有來源的輸入。

主要欄位：Mini/Mass/int_IMF/Gmag/G_BPmag/G_RPmag，或 initial_mass_msun/current_mass_msun/int_IMF/M_G/M_BP/M_RP。
程式對累積 int_IMF 作節點積分，不把等時線每列當等機率。
另一種明確命名的近似為 `--weight-mode power_law_initial_coordinate --alpha ...`；這不是自動測得的現今逃逸質量函數。

質量基準：

- `isochrone_current_stellar_mass`：提供範圍內的當前恆星質量，驗證 sum(current_mass * number_per_msun)=1。
- `initial_formed_mass`：初始形成質量等價的正規化，不能直接稱為今天逃離的恆星質量。

程式不自動補充未涵蓋的低質量星、殘骸、雙星或多族群。這些都要在人口來源與敏感度中說明。
外部已校正的加權表可透過 `Population(table,metadata).save(...)` 匯入；仍需相同欄位、質量正規化及來源。

輸出 population.csv＋population.json；包含模型、質量界限、每單位質量星數、mass_basis 與 hash。

## 4. 準備未做成員切選的 Gaia 母樣本

先生成可檢查的 ADQL：

```powershell
.\.venv\Scripts\python.exe selection_model.py parent-query --ra 229.02 --dec -0.121 --radius 1 --g-max 21 --out outputs/quality_parent/query.adql
```

若要從 ESA TAP 直接下載，使用新路徑：

```powershell
.\.venv\Scripts\python.exe selection_model.py parent-query --ra 229.02 --dec -0.121 --radius 1 --g-max 21 --download --max-records 250000 --out outputs/quality_parent/parent.csv.gz
```

此查詢只限制天區與 G，包含科學切選所需欄位、誤差相關係數與 G/BP/RP flux SNR。
不預先加 PM、視差、RUWE 或 CMD 切選。下載器檢查 VOTable QUERY_STATUS，遇到 OVERFLOW 或錯誤拒絕把截斷資料當完整母樣本。
大天區應分小 patch／採用外部可驗證的完整匯出，不建議直接下載密集區的大圓。
使用外部匯出時，metadata 必須保存實際 ADQL、來源說明及 parent_is_unfiltered=true。

球面半徑 1 度只示範取資料方法，不足以校正 5 度分析天區的所有位置；校正 grid 必須覆蓋真正使用範圍。

## 5. 校正缺測與聯合品質 map

grid.json 介面：

```json
{
  "center": [229.02, -0.121],
  "x_edges": [-1, 0, 1],
  "y_edges": [-1, 0, 1],
  "g_edges": [12, 16, 18, 19, 20, 21],
  "c_edges": [-0.5, 0.5, 1, 2]
}
```

x/y 是相對中心的 gnomonic 座標（度）；c 是 G-RP。實際格子要按樣本量與支持範圍選擇，以上不是通用最佳設定。

```powershell
.\.venv\Scripts\python.exe selection_model.py calibrate-quality --parent outputs/quality_parent/parent.csv.gz --parent-metadata outputs/quality_parent/parent.csv.json --grid grid.json --min-parent 20 --out outputs/quality_calibrated
```

拆成兩個條件：

1. P(科學所需欄位可得 | Gaia,x,y,G)。
2. P(原 RUWE／flux-excess 聯合品質合格 | 可得,x,y,G,G-RP)。

G 科學上限不在品質 map 重複施加；它會在觀測抽樣後明確執行。
map 使用 n/k 與 Beta posterior；同一校正格的抽樣會由所有粒子共享。
外推、n=0 或低於 min-parent 的格子會報錯，不以 1 或先驗均值當成可靠完整度。

限制：缺顏色來源的可得率先以 sky/G 建模；向尾巴人口轉移的適用性仍須獨立驗證。
本程式沒有把「通過 GC 成員選星／所有 Gaia」當作尾巴完整度。

## 6. response.json 設定

所有檔案路徑相對於 response.json，不是終端機目前目錄。

```json
{
  "schema_version": 2,
  "population": "population_calibrated",
  "quality": "quality_calibrated",
  "mass_basis": "isochrone_current_stellar_mass",
  "completeness": {
    "backend": "tcg_map",
    "mode": "multi",
    "path": "cache/gaia/allsky_uniq_10.fits",
    "source_description": "Official TCG DR3 reference map"
  },
  "errors": {
    "path": "quality_parent/parent.csv.gz",
    "metadata": "quality_parent/parent.csv.json",
    "covariance_policy": "required",
    "feature_scales": [0.5, 0.3, 5, 5],
    "max_feature_distance": 2,
    "neighbours": 20
  },
  "extinction": {"a_g": 0.1, "a_bp": 0.13, "a_rp": 0.07},
  "parallax_zero_point_mas": 0,
  "zero_point_description": "Explicit baseline; replace or test against a documented Gaia zero-point treatment",
  "observation_draws": 32,
  "calibration_draws": 16,
  "seed": 42,
  "cache_dir": "cache/response",
  "source_description": "Document actual calibrated inputs and assumptions",
  "synthetic_calibration": false
}
```

消光及零點數字只是語法範例，不是 Pal5 的校正值。
宣告 false 並不會讓不正確資料變成校正；任一輸入 metadata 標示 synthetic，結果會持續標示 synthetic。

其他完整度後端：

- `gaiaunlimited`：使用原套件 class；須能載入原生相依套件，失敗時不會默默換模型。
- `tabulated`：NPZ 提供 x_edges/y_edges/g_edges/probability，供外部校正與軟體測試。

## 7. 相關誤差與觀測響應

要求 parallax_pmra_corr、parallax_pmdec_corr、pmra_pmdec_corr；檢查係數與協方差矩陣半正定性。
明確選 `diagonal_approximation` 才允許對角近似。
光度誤差來自 flux SNR 或 phot_*_mean_mag_error，不默默填成零。

誤差 donor 條件於原始品質通過樣本，按 G、G-RP、位置匹配；記錄輔助誤差欄位缺測與 donor 支持範圍。
品質機率已包含的 RUWE／flux-excess，不會再生成一次並重複扣除。
光度誤差目前無跨波段相關，Gaussian 擾動明確截在 8 sigma，避免對未支持的任意暗星偷偷外推。
查詢完整度／品質時以無噪聲 apparent photometry 近似可觀測光度，須測試在邊界附近的偏差。

對每個粒子及人口格，以代表質量乘每單位質量星數，再乘入表、可得、品質與科學入選響應，得到真正的期望計數。
人口表在光度節點上積分；observation_draws 控制測量抽樣。
calibration_draws 控制共享品質 map 的 Beta 抽樣；用相同觀測亂數減少比較噪聲。

## 8. 一鍵執行與重用

```powershell
.\.venv\Scripts\python.exe run_pipeline.py --cluster Pal5 --legacy --mode weighted --reuse-run outputs/Pal5_1gyr --response-config outputs/response.json --out outputs/Pal5_weighted
```

亦可單獨重算觀測響應：

```powershell
.\.venv\Scripts\python.exe response_model.py --selection outputs/Pal5_weighted/selection --run outputs/Pal5_1gyr --config outputs/response.json --out outputs/response_variant
.\.venv\Scripts\python.exe fit_rate.py --input outputs/response_variant/calibrated_counts.npz --out outputs/variant_rate.json
```

輸出 response_maps.npz、calibrated_counts.npz、response.json、response.png。
包含每階段期望保留數、MC standard-error map、校正抽樣、mass_basis、input/code hash、版本與 cache key。
資料與程式版本相同時可重用 cache，且會檢查完整性。cache 與輸出必須是不同且互不包含的目錄。

## 9. 流失率與模型比較

固定動力學、釋放歷史形狀時：lambda_i = A*S_i + beta*B_i；rate=A*reference_rate。
不用每個 A 重跑動力學。N_i 為原始整數計數，S_i 是預期星數，不是 tracer 數或平滑背景扣除圖。

若有共享校正抽樣，對每個完整資料集的 likelihood 作平均，再取 log。
不逐格混用不同的校正實現。輸出的區間仍是 local asymptotic profile，在零邊界、低計數與多峰情況需另作檢查。
它自動納入的是提供的品質／可得率 map 抽樣；沒有憑空包含人口、銀河勢、背景錯配等未提供的系統誤差。

```powershell
.\.venv\Scripts\python.exe compare_models.py --templates outputs/modelA/calibrated_counts.npz outputs/modelB/calibrated_counts.npz --out outputs/model_comparison
```

程式核對觀測格子、背景、選星與質量定義一致，回報相對 likelihood 分數與參數。
delta_twice_nll 不是自動的偵測顯著性或 Bayes factor；掃描模型的試驗效應另需校正。

## 10. 時間變動流失率接口

CSV 欄位：lookback_myr、escaped_stellar_rate_msun_myr。必須覆蓋 [0,duration]，流失率非負。
它描述進入星流的恆星逃逸質量，不能直接拿包括氣體失質的總 dM/dt 代替。

在 response config 加：

```json
{
  "rate_history": "escape_history.csv",
  "rate_history_mass_basis": "isochrone_current_stellar_mass",
  "rate_history_source": "Actual physical model, version, component and units"
}
```

若重用尚無 release_time_sampling 標記的舊模擬，需明確加 assume_uniform_release_for_legacy=true。
對均勻釋放粒子，權重為 duration/N * rate(lookback_release_time)。回報的 reference_rate 是歷史在模擬期間的時間平均。
有限樣本的加權質量不強迫等於解析積分，避免引入不必要的正規化偏差。

這是固定動力學的再加權接口，沒有內建 EMACSS，也不會自動使星團 M(t)、半徑或釋放條件回饋演化。
顯著質量演化需另建一致的動力學模型；目前不能宣稱這項已完成。

軟體示範可加 `examples/create_v2_demo.py --variable-rate-demo`，產生明確標示 synthetic 的歷史。

## 11. 恢復與收斂檢查

```powershell
.\.venv\Scripts\python.exe validate_response.py --input outputs/response_variant/calibrated_counts.npz --amplitudes 0 0.5 1 2 --trials 20 --out outputs/response_recovery
.\.venv\Scripts\python.exe -m unittest -v test_model test_selection test_extent test_observation test_dynamics
```

validate_response 在提供的觀測響應下產生 Poisson 假資料、重新估計幅度與區間，回報偏差和覆蓋率。
這是參數模型下的恢復測試，不是獨立的品質／完整度校正，也不會重新建立 raw-catalogue CMD 或搜尋方向。
原 validate_tail.py 仍保留為固定背景、選後形態測試。

正式推論還需要人口、消光、背景模型、零點、空間系統誤差與不同星團的獨立驗證。
建議增加 tracer 數、observation_draws、calibration_draws 並更換 seed，確認結果穩定。

## 12. 本輪的積分效能改進

純 Python 預設改為 Cartesian DOP853，快取星團 Cartesian 軌跡 spline，直接評估原 MW forces 與解析 Plummer force。
僅適用目前 axisymmetric MW＋球對稱 Plummer 模型，沒有改變 Chen 釋放分布。
`stream_model.py --python-integrator galpy_reference` 保留原路徑作對照。

力評估已對照 galpy 的移動 Plummer 勢。相同 seed 的 1 Gyr／10 粒子案例，釋放時間完全相同，
天空位置差約 7e-6 度、自行差約 4e-7 mas/yr；原路徑約 183 秒，新路徑約 16 秒。
這是特定案例；run.json 和 manifest 會記錄實際耗時，不把小樣本速度直接當作所有模型的保證。

較大範例也已實跑：Pal5／1 Gyr／2,000 tracer 的模擬約 194 秒，包含選星、四組品質校正抽樣與擬合的整套流程約 214 秒。
這次使用真實 Gaia field 與官方完整度地圖，但人口／品質／誤差校正為明確標記的合成資料，因此只有軟體驗證意義。
本輪自動測試共 24 項；實際也測試了變動 rate 歷史、同觀測資料的模板比較，以及 cache 命中與完整性核對。
另外已從 ESA TAP 實際下載 0.02 度小天區的 13 筆 Gaia 母樣本，驗證欄位與 QUERY_STATUS；該小樣本不是正式品質校正。
觀測響應的參數恢復已實跑 amplitude=0/1 各 20 組。樣本量仍小，尤其零邊界的漸近區間不能假定具有名目覆蓋率；詳見輸出的二項區間及限制。
