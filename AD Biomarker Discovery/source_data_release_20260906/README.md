# SEA-AD source-data and statistical reproduction package

This package supplements, rather than replaces, the existing April discovery-run
archive. It supports the three representative CA1 reactive-astrocyte measurements,
the accumulated-panel calibration, the four corrected WSI baselines, the independent
GFAP comparison, and the aggregate external-cohort figures.

## Contents

- `inputs/run_*_train.csv` and `inputs/run_*_test.csv`: reconstructed, original-stack
  donor measurements and covariates for the five April searches. Machine-specific
  paths were replaced with slide filenames; numerical inputs are unchanged.
- `inputs/seaad_donor_split.csv`: the 35 discovery and 34 held-out donors.
- `inputs/baselines/`: the 34 held-out predictions per corrected baseline and compact
  fitting metadata. These are from the standardized rerun, not the historical
  CLIP-feature analysis that was labelled UNI.
- `inputs/seaad_gfap_heldout.csv`: the two CA1 GFAP regions and their mean, joined to
  the same 34 held-out donors. Source: SEA-AD June 30, 2026 quantitative neuropathology
  release, `Hip-MEC` sheet, `percent GFAP positive area`, CA1-L and CA1-R.
- `inputs/seaad_neurofm_core_readouts.csv`: the three replacement-stack measurements
  on the same held-out donors.
- `inputs/seaad_followup_donors.csv`: SEA-AD donor observations underlying its boxplots.
- `statistics/`: all 14 retained spatial entries, component-level correction,
  formula-specific bootstrap intervals, leave-one-donor-out summaries, paired
  baseline comparisons, panel calibration, and the plotted/recomputed curve audit.
- `illustrations/`: selected ED6d example contours and the corresponding subfield
  coordinates for H20.33.018. These are illustrative data, not another validation
  cohort. Whole-slide pixels and the full-slide ED6a cell overlay are not bundled here;
  use the original slide release and the repository's bulk reproducibility bundle.
- `scripts/`: figure reconstruction, optional statistical recomputation, and the
  numerical analysis scripts used for the archived inferential outputs.
- `MANIFEST.sha256`: checksums of all other files in this package and the accompanying
  SEA-AD panel CSVs, relative to the repository root.

## Measurement and statistical conventions

The outcome is `slope_zmem0`; covariates are age (`max_age_vis`), sex, Braak stage,
and CERAD. Each regression includes an intercept. More-negative signed correlations
indicate faster memory decline with greater reactivity. Direction-aligned plots
multiply the three tissue readouts' signed correlations by -1; baseline predictions
retain their natural direction.

At each discovery round, accumulated measurements are fitted separately within each
cohort. Multiple partial R is a within-cohort fitted joint association, not the
accuracy of unchanged discovery-trained panel weights on held-out donors. Missing
measurements are handled by complete-case analysis, not silently imputed. In run 2,
two hotspot-derived discovery columns have seven missing entries; its later-round
complete-case n is 28, while its discovery cohort still contains 35 donors.

The held-out panel null uses 100,000 shared-donor Freedman-Lane permutations (seed 11).
The round-10 CI uses 10,000 shared-donor basic-bootstrap resamples (seed 0). Individual
biomarker intervals are donor-percentile intervals. The 14-entry spatial family is
tested using maxT, followed by Holm correction across the bulk and spatial claims.
Fourteen entries represent thirteen unique formula definitions, including an
independent rediscovery; they are not fourteen independent claims.

Paired baseline differences use 10,000 donor bootstraps and 100,000 paired
randomizations, with Holm correction across four baselines. The historical analytic
`partial_r_p_two_sided` column computed with n-2 degrees of freedom is deliberately
not included. The inferential p-values here are the explicitly named permutation tests.

The canonical original-stack donor tables are an April-compatible reconstruction:
the historical evaluator did not retain the exact donor-level Test tables, and a
slide's classification state subsequently changed. They reproduce the current
original-stack component results, but must not be described as byte-identical to
every original evaluator input. Historical aggregate files remain unchanged.

See [the panel README](../../source_data/README_SEAAD.md) for the small but real
differences among the original archive, current plotting workbook, and reconstructed
trajectory. In particular, 0.054 is the historical round-10 run SD, not the
reconstructed curve's exact run SD of 0.0522902362.

## Reproduce the numerical figures

From a checkout of this repository, install `scripts/requirements.txt`, then run:

```bash
python "AD Biomarker Discovery/source_data_release_20260906/scripts/reproduce_figures.py" --output-dir /tmp/seaad_figures
```

This reproduces the numerical panels from the public CSVs. The follow-up plot shows
the published box summaries; restricted ADNI donor-level jitter points are omitted.

To recompute the held-out panel null, panel interval, family tests, and paired
baseline comparisons from public SEA-AD inputs:

```bash
python "AD Biomarker Discovery/source_data_release_20260906/scripts/recompute_statistics.py" --output-dir /tmp/seaad_statistics
```

The default settings reproduce the declared 100,000-permutation/10,000-bootstrap
analyses. `--quick` reduces resampling for a smoke test and does not reproduce the
reported p-values or intervals. Do not replace the archived results with quick-mode
outputs. Minor floating-point differences can occur across library versions.

## External-cohort source data

The ADNI files in `source_data/` contain published aggregate values only. Donor-level
ADNI inputs are not redistributed, consistent with the manuscript's data-use statement.
Use approved ADNI access for individual observations. The dementia figures were
exported from the archived plotting-input JSON; this release does not claim a fresh
independent reconstruction of the clinical dementia labels from ADNI source tables.

The optional `recompute_external.py` accepts an authorized local analysis CSV. It
expects the frozen readouts and precomputed outcomes; it does not perform classifier
training, feature selection, or a new discovery search. The follow-up groups are
nested sensitivity analyses, not independent cohorts.

ADNI-MEM slopes were estimated in the parent cohort using a linear mixed-effects
model; they must not be replaced with within-donor OLS slopes. Dementia status used
the latest valid global CDR >=1, with 0 and 0.5 classified as no dementia. The SEA-AD
dementia comparison used the cohort's binary cognitive-status field.

## Sources

- [SEA-AD data portal](https://brain-map.org/consortia/sea-ad/our-data)
- [Existing April search archive](../discovery_runs)
- [ADNI data access](https://adni.loni.usc.edu/)
- [Official UNI model](https://huggingface.co/MahmoodLab/UNI)
- [NeuroFM](https://arxiv.org/abs/2512.05993)

No ADNI participant records, API keys, absolute workstation paths, raw whole-slide
images, or large foundation-model embeddings are included in this source-data package.
