# Ten-cluster transfer preparation

The requested transfer replaces ESO 452-SC11, NGC 1851, M92 and M30 with
M62, NGC 6304, NGC 6355, NGC 6541 and M2. Existing observations are retained.

| Target | Canonical ID |
|---|---|
| M4 | NGC6121 |
| NGC 6139 | NGC6139 |
| NGC 6366 | NGC6366 |
| NGC 6397 | NGC6397 |
| NGC 6535 | NGC6535 |
| M62 | NGC6266 |
| NGC 6304 | NGC6304 |
| NGC 6355 | NGC6355 |
| NGC 6541 | NGC6541 |
| M2 | NGC7089 |

Preparation is separate from full numerical analysis. The transfer producer
copies completed fields and their metadata; it waits for the earlier field
download pass before starting additional missing-field downloads. Selection
uses the existing external selector. Its required source files and catalogue
are copied into the transfer directory. Do not publish those private inputs
with the public code repository.

The model preparer uses an isolated source snapshot, a genuine uncut Gaia DR3
calibration parent, an additional bright-star census, the official GaiaUnlimited
order-7 completeness map and the supplied PARSEC isochrone. The downloaded
sky/PM-window science field is not used as an uncut calibration parent.
Unsupported quality-map queries remain errors; unused sparse colour cells do
not by themselves establish failure of the actual model response. This is a
support rule, not a validation of field-to-tail population transfer.

All 15 population nodes (three escape-slope offsets and five measured-MF
quadrature nodes each) must have positive weights and unit represented living
stellar mass. A small real-input integration and response must pass before
export. A separate subprocess then verifies the exported files and tests the
rebound paths against the same small simulation. A package is marked ready only
after that check. Full-duration response support, sampling convergence,
background checks and scientific recovery tests remain necessary.

New packages use `analysis_kind=model_comparison`: background and measured-MF
uncertainty are integrated, and the output includes relative likelihood grids
and profiles. No posterior distribution over rate and lookback time is
computed. Existing legacy packages retain their original analysis mode.
Confidence levels and a validated physical total mass-loss rate are not
certified by completing a numerical grid.

The single-isochrone baseline requires particular caution for M62's helium
populations ([Milone 2015](https://doi.org/10.1093/mnras/stu2198)) and M2's
minority metal-rich population ([Lardo et al. 2016](https://doi.org/10.1093/mnras/stv2802)).
These are conditional approximations, not models of all population components.
The escape-MF offsets also remain assumptions rather than measurements of the
tail MF. The population's mass basis excludes unmodelled remnants, winds and
stars below the supplied isochrone's mass limit.

For the current preparation, inspect these files under
`outputs/other_pc_ten_20261101_v1`:

- `transfer_plan.json`: fixed requested targets and analysis mode.
- `status.json`: observation preparation, including explicit failures.
- `model_preparation_status.json`: calibration, tests and package readiness.
- `packages/*_portable_v1.zip`: packages that have reached the export stage.
- `package_validation/<ID>/validation.json`: required relocation-check result.

The existence of a ZIP alone is not sufficient: export precedes the relocation
test. Only distribute a package whose model status is
`conditional_model_package_ready`. Each package includes its own `RUN_GUIDE.md`,
checksums and dependency versions. Raw calibration-parent duplicates, tile
caches, environments and old numerical results are excluded from the export;
the empirical error-donor catalogue required at runtime is included.

The observation producer waits for at least 4 GiB of available system RAM
before selecting a field; model preparation waits for at least 8 GiB before
each substantial stage. These gates avoid starting new heavy stages during
the existing local numerical run; they are not hard memory limits.

After stopping any earlier writers, restart the same preparation commands in
the repository's existing environment to reuse compatible completed stages:

```powershell
.\.venv\Scripts\python.exe -u prepare_transfer_batch.py --root outputs/other_pc_ten_20261101_v1 --previous-transfer outputs/other_pc_nine_20261101_v1
.\.venv\Scripts\python.exe -u prepare_transfer_models.py --root outputs/other_pc_ten_20261101_v1
```

Run the model preparer in a second terminal if preparation should overlap the
observation producer. Do not start two copies of either producer against the
same output directory. Interrupted incomplete scientific outputs are preserved;
review failures before choosing a fresh repair version.
