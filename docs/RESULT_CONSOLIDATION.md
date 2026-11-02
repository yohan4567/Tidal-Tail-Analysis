# Consolidating completed results

`consolidate_completed_results.py` builds a new result snapshot without moving
source results, editing frozen inputs, starting dynamics, or changing active jobs.
It reads the local comparison folders and discovers every numeric `NGC` folder
under the adjacent `results` directory, including newly returned clusters.

```powershell
.\.venv\Scripts\python.exe consolidate_completed_results.py --out outputs/consolidated_results_NEW_VERSION --doc "C:/path/to/Claude資料表.docx"
```

Use a fresh output directory. The supplied document contains the 29-cluster
physical-parameter and morphology tables; its rounded physical values must agree
with the project's frozen catalogue. Morphology labels are recorded as supplied,
not independently validated by this exporter.

The Python command produces CSV/JSON tables, canonical comparison grids and
profiles, verified source copies, and a manifest. The formatted Excel workbook in
the November 2 snapshot was authored separately; this command does not require
or install a spreadsheet runtime and does not create that workbook.

Only clusters with all three offsets enter the completed-result table. Score-only
returns are preserved under `pending_comparison`, outside the integrated
likelihood summary. Duplicate local/returned versions cause an error instead of
silently choosing a version.

Checks include source-copy SHA256 agreement, raw JSON/CSV equality, complete
16 x 20 x 5 nodes, offset identity, recomputed q and profiles, and reconstruction
of the integrated surface from the background-integrated nodes and population
weights. For returned runs the original runtime population-prior file may be
absent: compatible local package weights reproduce the surface, but this does
not establish the missing runtime file's hash or complete execution provenance.

The MLR estimate remains a sampled grid point. The q <= 1 span is a piecewise
linear diagnostic, not a calibrated confidence interval. Truncated directions
remain unavailable. MLR/M uses the adopted present-day catalogue mass and does
not propagate its error; the represented population mass and catalogue total
mass are not identical physical scopes.
