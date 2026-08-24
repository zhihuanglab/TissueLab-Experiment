# TissueLab-Experiment

## Project Overview

This repository contains experimental data and video demonstrations for the TissueLab project, showcasing how to visualize outputs and experimental results across multiple medical image analysis tasks. The project covers cell counting, cell proportion analysis, depth of invasion measurement, fatty liver detection, hypertrophy detection, intracranial hemorrhage identification, kidney glomerulus counting, lymph node counting, metastasis classification, and X-ray image analysis.

## 🖥️ YouTube Demonstration Video
- [TissueLab Demonstrations](https://www.youtube.com/watch?v=rssWT4Mehqw) - A demonstration video showcasing TissueLab experimental results and visualizations

## 🔬 Reproducibility Bundles

Reproducibility bundles for the AD Biomarker Discovery, Metastasis classification and Tubule Score discovery experiments are hosted on Google Drive. The bulk intermediate outputs of the external validation experiments are hosted there as well, since they are too large for this repository:

- [TissueLab Reproducibility Bundles (Google Drive)](https://drive.google.com/drive/folders/19YExomssX5Pz3D7Jg1yJCUU9uhMQc7FM?usp=sharing)

## Video Demonstrations
**Recordings/** - Contains demonstration videos showing how to visualize outputs and reproduce experiments for better replication:

- `depth_of_invasion.mov` - Depth of invasion measurement demonstration
- `fatty_liver.mp4` - Fatty liver detection demonstration  
- `hemorrhage.mp4` - Intracranial hemorrhage identification demonstration
- `integrate_any_model_as_tasknode.mp4` - Model integration demonstration
- `lymph_node_counting.mp4` - Lymph node counting demonstration
- `refine_classifier_and_apply_all.mp4` - Classifier refinement and application demonstration

## Experimental Data Structure

### Dataset Folders

- **Cell Counting/** - Neoplastic cell quantification. A pathologist taught the
  system to recognise neoplastic cells on Visium HD colon tissue; the resulting
  classifier was then frozen and carried unchanged to a Xenium colon section, where
  the reference is defined transcriptomically rather than morphologically.
  - Q1.7round1.csv to Q1.7round5.csv (5 rounds) - the Visium HD source platform.
  - Q1.7_Xenium_round1.csv to Q1.7_Xenium_round5.csv (5 rounds) - the unseen Xenium
    section. `ground_truth` is 261,905 neoplastic cells; `TLAgent` answered 253,870
    in every run, `TLAgent-removing knowledge` is the same pipeline without the
    expert annotations.
  - `Xenium_colon_per_cell.parquet` - the transcriptomic reference and the
    classifier score for each of the 446,933 cells on the section: the cell id, its
    centroid, its transcriptomic cell type, whether that type is neoplastic, and the
    predicted probability. Scoring `prob_positive` against `gt_positive` gives the
    reported cell-level area under the ROC curve.

- **Cell Proportion/** - Prostate epithelial cell proportion. The same recognition
  module taught on Visium HD prostate tissue, frozen and applied to an unseen Xenium
  prostate section. The quantity asked for differs between platforms because the two
  references support different quantities: the tumour-to-duct ratio on Visium HD and
  the proportion of normal duct cells among all cells on Xenium. What is carried
  across unchanged is the recognition module, not the downstream computation.
  - Q1.6.1round1.csv to Q1.6.1round5.csv (5 rounds) - the Visium HD source platform.
  - Q1.6.1_Xenium_round1.csv to Q1.6.1_Xenium_round5.csv (5 rounds) - the unseen
    Xenium section. `ground_truth` is 47.15 per cent and `TLAgent` answered 42.09
    per cent in every run.
  - `Q1.6.1_Xenium_extval_normalduct.csv` - the answer averaged over the five runs
    with its absolute error against the transcriptomic reference.
  - `Xenium_prostate_per_cell.parquet` - the transcriptomic reference and the
    classifier score for each of the 156,021 typed cells, in the same layout as the
    colon file. Scoring `prob_positive` against `gt_positive` gives the reported
    cell-level area under the ROC curve.

- **Depth of Invation/** - Depth of invasion measurement task data (internal LNCO2 cohort)
  - Q1.1round1.csv to Q1.1round5.csv (5 rounds of experimental data). 107 rows;
    the 2 slides carrying `ground_truth = N/A` are those on which the depth could
    not be determined with confidence, leaving 105 evaluable cases from 43 patients.
  - GT_LNCO2_doi_annotation_lines.csv - reference standard. Our collaborating
    pathologist measured the invasion depth directly on each whole-slide image by
    drawing a line from the surface to the deepest point of invasion. Endpoints in
    slide pixel coordinates, the slide's microns per pixel, and the resulting
    length in millimetres. 105 rows. `Length_mm` is the `ground_truth` column of
    the round files; join on `ImageName` = `slide_name`.
  - LNCO2_splits/ - slide-to-patient mapping and split membership, needed to
    reproduce the patient-level bootstrap intervals since several patients
    contribute many slides. `lnco2_split_slides.csv` gives `patient_id` and the
    `in_TRAIN` / `in_Q1.1_DOI` / `in_Q1.2_LNcount` / `in_Q1.3_Metastasis` flags for
    476 slides; `lnco2_split_patients.csv` gives per-patient counts for 50 patients.
    The evaluation sets share no slide and no patient with the set used for expert
    feedback. Covers Q1.1, Q1.2 and Q1.3, which all draw on the same cohort.

- **Depth of Invation/External Validation TCGA-COAD/** - the workflow co-evolved on
  LNCO2 was frozen and applied once to TCGA-COAD, with no further co-evolution
  rounds, no new annotation and no fine-tuning. 82 whole-slide images from 77
  patients contributed by seven independent tissue source sites (TSS codes 3L, 4N,
  AZ, CK, G4, NH, QG), which differ from LNCO2 in institution, country, scanner and
  staining protocol.
  - Q1.1round1.csv to Q1.1round5.csv (5 rounds of experimental data). 82 rows;
    the pathologist could determine the reference on 65 of them, from 63 patients
    and spanning all seven sites, and the remaining 17 carry `ground_truth = N/A`.
    All reported statistics are computed on the 65 evaluable cases. The annotation
    was completed in full before TissueLab or any baseline was run on this cohort.
  - TCGA_COAD_cohort.csv - one row per slide, all 82, with `patient_id`, `tss_code`,
    `evaluable`, `invasion_depth_mm` and `exclusion_reason`.
  - annotation_lines.csv - reference standard, same format as the internal cohort.
    66 rows for 65 slides: `TCGA-CK-4947-01Z-00-DX1` carries two lines drawn 77
    seconds apart, the later superseding the earlier. Where a slide carries more
    than one line, take the one with the greatest `datetime`.
  - Methods differ substantially in how often they return a valid numeric answer
    on this cohort, so error metrics computed over different subsets are not
    directly comparable; count the non-`N/A` entries per method before comparing.

- **Fatty Liver/** - Fatty liver detection task data
  - Q2.1round1.csv to Q2.1round5.csv (5 rounds of experimental data)

- **Hypertrophy/** - Hypertrophy detection task data
  - Q2.6round1.csv to Q2.6round5.csv (5 rounds of experimental data)

- **Intracranial Hemorrhage/** - Intracranial hemorrhage identification task data
  - Q2.4.1round1.csv to Q2.4.1round5.csv (5 rounds of experimental data)
  - Q2.4.2round1.csv to Q2.4.2round5.csv (5 rounds of experimental data)

- **Kidney Glomerulus Counting/** - Kidney glomerulus counting task data
  - Q3.2round1.csv to Q3.2round5.csv (5 rounds of experimental data)

- **Lymph Node Counting/** - Lymph node counting task data (LNCO2 cohort)
  - Q1.2round1.csv to Q1.2round5.csv (5 rounds of experimental data). 321 adjacent
    lymph node slides from 22 patients, between 6 and 28 slides per patient.
    `ground_truth` is the number of metastatic lymph nodes on the slide, an integer
    from 0 to 5. The cohort is imbalanced, 280 of the 321 slides carry none, which
    is why the weighted F1-score is reported alongside overall accuracy. This
    reference comes with the LNCO2 collection rather than being annotated for this
    study.
  - `GT_LNCO2_lymph_node_annotations.csv` - the region-level annotations the counts
    are read from, one row per slide. `n_positive_LN` is the `ground_truth` column
    of the round files and `n_total_LN` is how many lymph nodes the slide contains;
    `coordinates` holds the annotated regions as normalised polygon vertices, with
    `name` giving the region type for each. `case_name` is the patient. Join on
    `slide_name`.
  - Slide-to-patient mapping is in `Depth of Invation/LNCO2_splits/`
    (`in_Q1.2_LNcount` column), needed to reproduce the patient-level bootstrap
    intervals. These slides share no patient with the set used for expert feedback.

- **Metastasis Classification/** - Metastasis classification task data (LNCO2 cohort)
  - Q1.3round1.csv to Q1.3round5.csv (5 rounds of experimental data). The same 321
    lymph node slides from 22 patients as the counting task. `ground_truth` is the
    metastasis category: macrometastasis (> 2.0 mm), micrometastasis (>= 0.2 mm),
    isolated tumour cells, or negative. The `TLAgent-R1`, `TLAgent-R2`,
    `TLAgent-R1+2` and `TLAgent-R1+2+Reflector` columns are the system after the
    corrective skills distilled from its own errors have been applied, singly and
    in combination.
  - `GT_LNCO2_metastasis_annotations.csv` - the reference standard. The dataset
    marks which sections contain tumour-positive areas but does not give the deposit
    size, so our collaborating pathologist measured each deposit on the whole-slide
    image and assigned the category. One row per measured deposit, 70 rows over the
    37 slides that carry one, with endpoints in slide pixel coordinates, the slide's
    microns per pixel, and the length in millimetres. `slide_category` is the
    `ground_truth` column of the round files; join on `ImageName` = `slide_name`.
    The category follows from the largest deposit on the slide, above 2.0 mm for
    macrometastasis and at or above 0.2 mm for micrometastasis.
  - `Q1.3_planner_llm_variation.csv` - planner backbone sensitivity. Four backbones
    each generated 50 independent plans for this task with no other component
    changed, 200 rows. Columns record the plan, its step count, and whether it
    contains each structural property; `error` is non-empty only if no executable
    plan was produced, which did not occur.
  - `Q1.3_reflector_llm_variation_round1.csv` to `..._round5.csv` - reflector
    backbone sensitivity. Per-slide predictions on the same 321 slides when the
    VLM-assisted reflector is driven by each of four backbones, over five repeated
    runs.
  - Slide-to-patient mapping is in `Depth of Invation/LNCO2_splits/`
    (`in_Q1.3_Metastasis` column).

- **Tubule Score/** - Nottingham tubule formation grading on TCGA-BRCA, an ordinal
  1 to 3 scale over 30 test cases. Five independent agents each searched for an
  analytical workflow over 20 rounds.
  - Q1.5round1.csv to Q1.5round5.csv - the five searches at their final round,
    alongside four training-based baselines and three human-designed workflows.
  - `data/agent_1/` to `data/agent_5/` - one directory per search.
    `agent_N_summary.csv` gives the metrics at each of the 20 rounds, and
    `eval/test/round{r}_predictions.csv` the per-case predictions that produced them.
  - `data/candidate_configs/` - every configuration the search generated at each
    round, selected or not, with its leave-one-out correlation, its in-sample
    correlation, the gap between them and the penalised score that ranked it. The
    `round{r}_is_predictions.csv` files hold the leave-one-out predictions behind
    those scores. This is what makes the selection behaviour replayable rather than
    something to be taken on trust.
  - `data/ablation/` - the same search restricted to a fixed simple classifier
    (`suboptimal/`) against the unrestricted search (`optimal/`), five agents each.
  - `data/baselines/` - the four training-based baselines.
  - `data/human_workflows/` - the three human-designed workflows.
  - `External Validation histAI/` - the workflows discovered on TCGA-BRCA applied
    without further rounds to 75 cases from HISTAI, an independent collection
    prepared and scanned elsewhere. No HISTAI label entered the search, and the
    external cohort was scored only after each round had already been selected on
    internal data.
    - `histai_external_predictions.csv` - the reference grade and each of the five
      agents' continuous prediction for all 75 HISTAI cases. Converting these to
      per-class scores with a Gaussian softmax over the class centres, as the
      Methods describe, gives the reported external macro-AUC.
    - `tcga_test_predictions.csv` - the same for the 30 internal test cases.
    - `histai_external_round_curve.csv` - the external correlation at each of the
      20 co-evolution rounds, per agent and pooled.
    - `three_split_summary.csv` reports each agent in sample on its training set, on
      leave-one-out, on the TCGA test set and on HISTAI; the other files compare the
      agents against the baselines on both cohorts.

- **X-Ray/** - X-ray image analysis task data
  - Contains detection data for multiple pathological types:
    - Atelectasis
    - Cardiomegaly
    - Consolidation
    - Edema
    - Effusion
    - Emphysema
    - Fibrosis
    - Hernia
    - Infiltration
    - Mass
    - Nodule
    - Pleural_Thickening
    - Pneumonia
    - Pneumothorax
  - Each pathological type contains 5 rounds of experimental data
  - Also includes co-evolving related data files

## Data Format

CSV files contain the following columns:
- `slide_name` - Slide name
- `ground_truth` - Ground truth labels/values
- `TLAgent` - TLAgent model prediction results
- `TLAgent-removing knowledge` - the identical system evaluated in its state
  before co-evolution, that is with the expert feedback removed and everything
  else unchanged. Present where the paper reports a before-and-after comparison.
- `Quilt-LLaVA` - Quilt-LLaVA model prediction results
- `LLaVA-Med` - LLaVA-Med model prediction results
- `MedGemma` - MedGemma model prediction results
- `GPT-4o-agent` - GPT-4o-agent model prediction results
- `GPT-4o-vision` - GPT-4o-vision model prediction results
- `GPT-5-vision` - GPT-5-vision model prediction results

## Notes

- Some models may show "N/A" for certain tasks, indicating the model was unable to complete the task
