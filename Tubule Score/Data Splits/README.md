# Data Splits

Two stand-alone scripts that regenerate every case split used in the Tubule
Score experiments. Each reads only the pool file next to it, uses the Python
standard library, and writes its split CSVs (`case_id`, `tubule_score`) into
this folder. Both use seed 42.

    python split_tcga_brca.py   # -> tcga_brca_train.csv, tcga_brca_test.csv
    python split_histai.py      # -> histai_batch{1,2,3}.csv, histai_external75.csv

| Output | n | Per score | Used for |
|---|---|---|---|
| `tcga_brca_train.csv` | 30 | 10 / 10 / 10 | agent search, `Training LOO TCGA-BRCA/` |
| `tcga_brca_test.csv` | 30 | 10 / 10 / 10 | `Internal Validation TCGA-BRCA/`, `TCGA-BRCA/Q1.5round*.csv` |
| `histai_batch1.csv` | 30 | 10 / 10 / 10 | first external accrual batch |
| `histai_batch2.csv` | 30 | 10 / 10 / 10 | second batch, disjoint from batch 1 |
| `histai_batch3.csv` | 15 | 5 / 5 / 5 | third batch, disjoint from 1 and 2 (score 1 exhausted) |
| `histai_external75.csv` | 75 | 25 / 25 / 25 | `External Validation histAI/` |

Inputs: `tcga_brca_pool.csv` (542 TCGA-BRCA cases with a report-stated tubule
score, with the supporting excerpt in `note`) and `histai_breast_pool.csv`
(858 HISTAI-breast cases with a tubular-differentiation score in the metadata).

Run without arguments, each script also checks its regenerated split against
the per-case prediction files in the sibling folders and fails on any mismatch.
The exact sampling recipe is documented at the top of each script.
