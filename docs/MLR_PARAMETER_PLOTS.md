# Delta-zero MLR and cluster-parameter plots

`plot_mlr_parameters.py` reads completed delta=0 integrated-likelihood comparisons
and the verified legacy integrated surfaces for NGC 288, Pal 5 and M107. It does
not compute a posterior or combine score-only maxima with integrated results.
Input score hashes, profiles, offsets and background integration settings are
checked. Clusters without a matching comparison are recorded as excluded.

The supplied DOCX contains tail-shape labels and rounded physical parameters.
All 29 physical rows are checked against the project's saved Baumgardt catalogue
at the DOCX precision. Plot coordinates retain the catalogue precision;
eccentricity is `(R_apo-R_peri)/(R_apo+R_peri)` from the catalogue orbit model.
Shape labels are reproduced from the supplied table, not independently validated.

Vertical spans show piecewise-linear `q_profile <= 1` regions. They are diagnostic
display intervals, not standard deviations, credible intervals or calibrated
confidence intervals. Grid-truncated directions are null in the numeric error
columns and uncapped in the plots. Empty best-fit markers flag grid-edge maxima.
The two-decimal CSV is a display export; it does not establish 0.01 precision.

```powershell
.\.venv\Scripts\python.exe plot_mlr_parameters.py --doc "path\to\reference.docx" --out outputs/mlr_parameters --fit-linear
```

Without `--fit-linear`, the mass panel uses a log display axis. With it, all three
panels fit `y = a*x + b` by unweighted OLS. The mass coordinate is linear
`M/(10^5 Msun)`, not log mass. The q spans are not used as weights. Fits include
grid-edge best points and are descriptive; no coefficient intervals, detection
significance or causal interpretation is supplied. Coefficients are checked
against the independent centred OLS formula and R-squared against squared
Pearson correlation. The saved data and source hashes accompany PNG/PDF/SVG
figures, the fit table, a combined figure and interpretation notes.
