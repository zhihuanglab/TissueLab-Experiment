# TissueLab-Experiment

## Project Overview

This repository contains experimental data and video demonstrations for the TissueLab project, showcasing how to visualize outputs and experimental results across multiple medical image analysis tasks. The project covers cell counting, cell proportion analysis, depth of invasion measurement, fatty liver detection, hypertrophy detection, intracranial hemorrhage identification, kidney glomerulus counting, lymph node counting, metastasis classification, and X-ray image analysis.

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

- **Cell Counting/** - Cell counting task data
  - Q1.7round1.csv to Q1.7round5.csv (5 rounds of experimental data)

- **Cell Proportion/** - Cell proportion analysis task data
  - Q1.6.1round1.csv to Q1.6.1round5.csv (5 rounds of experimental data)

- **Depth of Invation/** - Depth of invasion measurement task data
  - Q1.1round1.csv to Q1.1round5.csv (5 rounds of experimental data)

- **Fatty Liver/** - Fatty liver detection task data
  - Q2.1round1.csv to Q2.1round5.csv (5 rounds of experimental data)

- **Hypertrophy/** - Hypertrophy detection task data
  - Q2.6round1.csv to Q2.6round5.csv (5 rounds of experimental data)

- **Intracranial Hemorrhage/** - Intracranial hemorrhage identification task data
  - Q2.4.1round1.csv to Q2.4.1round5.csv (5 rounds of experimental data)
  - Q2.4.2round1.csv to Q2.4.2round5.csv (5 rounds of experimental data)

- **Kidney Glomerulus Counting/** - Kidney glomerulus counting task data
  - Q3.2round1.csv to Q3.2round5.csv (5 rounds of experimental data)

- **Lymph Node Counting/** - Lymph node counting task data
  - Q1.2round1.csv to Q1.2round5.csv (5 rounds of experimental data)

- **Metastasis Classification/** - Metastasis classification task data
  - Q1.3round1.csv to Q1.3round5.csv (5 rounds of experimental data)

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
- `Quilt-LLaVA` - Quilt-LLaVA model prediction results
- `LLaVA-Med` - LLaVA-Med model prediction results
- `MedGemma` - MedGemma model prediction results
- `GPT-4o-agent` - GPT-4o-agent model prediction results
- `GPT-4o-vision` - GPT-4o-vision model prediction results
- `GPT-5-vision` - GPT-5-vision model prediction results

## Notes

- Some models may show "N/A" for certain tasks, indicating the model was unable to complete the task
