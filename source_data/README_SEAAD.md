# SEA-AD figure source data

These files contain the numerical inputs to the SEA-AD discovery, measurement-stack
comparison, and external validation panels. They follow the repository's panel-level
CSV naming convention. The companion [analysis package](../AD%20Biomarker%20Discovery/source_data_release_20260906/README.md)
contains donor-level SEA-AD inputs, statistical outputs, and reproduction scripts.

| File | Panel | Contents |
|---|---|---|
| `extfig6b_seaad_panel_trajectory.csv` | ED6b | Current plotted discovery/test curves, run SD, held-out permutation null, and round-10 basic-bootstrap CI |
| `extfig6c_seaad_biomarker_comparison.csv` | ED6c | Four corrected WSI baselines, independently quantified GFAP, and three frozen tissue readouts |
| `extfig7k_seaad_adni_dementia.csv` | ED7k | Aggregate dementia-status associations and subsampling SD |
| `extfig8a_classifier_training_labels.csv` | Supplementary 4a | Cell-class label counts and exact class colors for both classifier versions |
| `extfig8b_seaad_tool_transfer.csv` | Supplementary 4b | Original/replacement-stack associations and subsampling SD |
| `extfig9a_cognitive_followup.csv` | Supplementary 5a | Aggregate boxplot quartiles and whiskers; no participant-level ADNI records |
| `extfig9b_adni_followup_associations.csv` | Supplementary 5b | Full ADNI, >=5-year and >=7-year groups, and rebuilt-stack SEA-AD reference |

## Column definitions

- `value` is the plotted statistic, not necessarily a signed association. Read `metric`.
- `signed_partial_r` preserves the original direction when `value` is direction-aligned.
- `sd_runs` is between-search SD. It is blank for once-fitted/fixed measurements.
- `sd_subsample` is SD over repeated donor subsets. It is not a standard error or a CI.
- `subsample_sd_ddof` records the implemented SD convention where verified. The archived
  external dementia aggregate did not retain that setting, so it is not invented here.
- `ci95_*` columns contain confidence limits only. A blank means not supplied, not zero.
- `n_cases` and `n_patients` describe donor cohort size. The ED6b per-run input table
  additionally records complete-case counts; the discovery runs do not all use 35
  complete cases at every round.
- `recomputed_value`, `recomputed_sd_subsample`, and `reconstruction_*` provide
  recalculated quantities alongside rounded or separately archived plotting inputs.
- `permutation_p_upper` is the upper-tail, pointwise fitted-panel permutation p-value;
  no corresponding discovery-cohort p-value is supplied.

## Important version distinction

ED6b's current plotted values come from `panel_trajectory_figure_data.xlsx`, sheet
`PerRunCurves`. The held-out null and CI come from the complete-case reconstructed
donor measurements. These are distinguished in the CSV rather than silently
replacing historical data. The older `aggregate_panel_train_test_5run_*.csv` files
in the experiment folder are historical archives, not the current figure inputs.

The reconstructed held-out round-10 R is 0.6890088234 with between-run SD
0.0522902362. The archived original trajectory has R 0.6894627344 and SD
0.0539164888. The current plotting workbook gives R 0.6890 and SD 0.0523015774
from four-decimal per-run inputs. Thus the manuscript's historical SD of 0.054
must not be described as the exact SD of the reconstructed curve.

The current discovery plotting workbook ends at 0.6790, whereas complete-case
recomputation from the supplied original-stack donor tables gives 0.6798706985.
See the companion `statistics/trajectory_version_comparison.csv`. The files
preserve this discrepancy for review instead of modifying values to match text.

## Access and privacy

Only aggregate ADNI results are included. The exact donor-level jitter points in
Supplementary 5a cannot be reconstructed from the public aggregate CSV. Authorized
ADNI users can supply their own approved participant-level inputs to the companion
recomputation script. Participant identifiers, clinical records, slide paths, and
participant-level ADNI biomarker values are not distributed in this package.
