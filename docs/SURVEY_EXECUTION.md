# Remaining-cluster queue

Run `python run_remaining_survey.py --workers 10 --response-workers 3 --prepare-workers 2`.
The default sample excludes completed Pal 5. Up to two preparation subprocesses
run while a single cluster grid uses the numerical worker budget. This prevents
multiplying ten dynamics workers by the number of clusters.

Per cluster, the queue uses its frozen observed selection, catalogue phase space,
mass/radius, literature age, its own official PARSEC request and unfiltered Gaia
calibration parent. Existing six-cluster Gaia parents are reused after the
calibration's checksum validation. Other supported clusters fetch their own
parents via ESA/AIP with resumable tiles. No Pal 5 population/error calibration
is transferred. Once ready, each cluster runs the same 16-rate, 20-time, three-seed
conditional grid with shared observation groups. Failures are recorded and the
queue continues.

`configs/survey_populations_v1.json` records adopted age sources and approximation
policies. These are conditional single-population baselines, not fitted CMDs.
Scaled-solar model [M/H] is numerically set to catalogue [Fe/H]; alpha enhancement,
helium, multiple populations and differential extinction are not calibrated by
that choice. Complex populations and missing documented ages are held for review.
The f=1 conversion from represented living-star escape to dynamical loss is an
explicit assumption. A completed grid is not a validated absolute escape rate.

Age references: VandenBerg et al. (2013), Table 2,
https://arxiv.org/abs/1308.2257 ; NGC2298: Forbes & Bridges (2010), age table,
https://academic.oup.com/mnras/article/404/3/1203/1049211 .

Progress artifacts, all under `outputs/cluster_survey_29_v1`:

- `PROGRESS.md` and `progress.json`: actual dynamics/templates plus queue stage.
- `remaining_survey_status.json`: preparation, ready, running, completed, blocked,
  or failed stages; detailed error/reason rather than inferred success.
- `<cluster>/remaining_job.json`: per-cluster state and adopted assumptions.
- `<cluster>/production_shared.log`: numerical grid launcher log; individual
  integration logs/checkpoints remain inside its grid directory.

Rerunning the command uses valid completed inputs and resumes the source-frozen
grid. Do not run two queue writers against the same root. Interrupted/inconsistent
calibration bundles must be inspected before retry; they are never overwritten.
Missing real observations are blocked instead of replaced with synthetic inputs.

## Model-grid boundary check

During this queue's first preparation, requests for M92 [M/H]=-2.31 and M30
[M/H]=-2.27 returned MH=-2.19174. Population construction rejected the mismatch.
The downloader now also checks returned age and MH before accepting its metadata.
Those two baseline preparations require a justified composition/model choice;
they must not be relabelled as the requested metal-poor models or passed by
loosening the numerical tolerance. Other cluster jobs continue independently.
