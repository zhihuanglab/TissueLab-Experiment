# Source Data workbooks

One workbook per figure, one sheet per figure panel, as the Nature manuscript checklist requires. Each sheet is the corresponding per-panel CSV in the parent folder, unchanged. Regenerate with `python3 build_workbooks.py` from `source_data/`.

| Workbook | Sheet | Source CSV | Rows |
|---|---|---|---|
| Source_Data_ED_Fig2.xlsx | ED Fig. 2c | extfig2c_ich_accuracy.csv | 10 |
| Source_Data_ED_Fig2.xlsx | ED Fig. 2d | extfig2d_ich_MCC.csv | 10 |
| Source_Data_ED_Fig2.xlsx | ED Fig. 2e | extfig2e_ivh_accuracy.csv | 10 |
| Source_Data_ED_Fig2.xlsx | ED Fig. 2f | extfig2f_ivh_MCC.csv | 10 |
| Source_Data_ED_Fig3.xlsx | ED Fig. 3b | extfig3b_xray_accuracy.csv | 84 |
| Source_Data_ED_Fig3.xlsx | ED Fig. 3c | extfig3c_xray_auc_before_after.csv | 28 |
| Source_Data_ED_Fig4.xlsx | ED Fig. 4a | extfig4a_tubule_run1_pearson_trajectory.csv | 40 |
| Source_Data_ED_Fig4.xlsx | ED Fig. 4b | extfig4b_tubule_comparison_pearson.csv | 8 |
| Source_Data_ED_Fig4.xlsx | ED Fig. 4c | extfig4c_tubule_baseline_training_dynamics.csv | 1204 |
| Source_Data_ED_Fig4.xlsx | ED Fig. 4d | extfig4d_tubule_classifier_ablation_loo_auc.csv | 40 |
| Source_Data_ED_Fig4.xlsx | ED Fig. 4e | extfig4e_tubule_classifier_ablation_test_auc.csv | 40 |
| Source_Data_ED_Fig6.xlsx | ED Fig. 6b | extfig6b_seaad_panel_trajectory.csv | 20 |
| Source_Data_ED_Fig6.xlsx | ED Fig. 6c | extfig6c_seaad_biomarker_comparison.csv | 8 |
| Source_Data_ED_Fig7.xlsx | ED Fig. 7a | extfig7ab_cellcounting_xenium.csv | 17 |
| Source_Data_ED_Fig7.xlsx | ED Fig. 7b | extfig7ab_cellcounting_xenium.csv | 4 |
| Source_Data_ED_Fig7.xlsx | ED Fig. 7d | extfig7de_cellproportion_xenium.csv | 18 |
| Source_Data_ED_Fig7.xlsx | ED Fig. 7e | extfig7de_cellproportion_xenium.csv | 4 |
| Source_Data_ED_Fig7.xlsx | ED Fig. 7g,h | extfig7ghi_doi_external.csv | 30 |
| Source_Data_ED_Fig7.xlsx | ED Fig. 7i | extfig7ghi_doi_external.csv | 4 |
| Source_Data_ED_Fig7.xlsx | ED Fig. 7j | extfig7j_tubule_transfer_auc.csv | 40 |
| Source_Data_ED_Fig7.xlsx | ED Fig. 7k | extfig7k_seaad_adni_dementia.csv | 6 |
| Source_Data_ED_Fig8.xlsx | ED Fig. 8a | extfig8a_classifier_training_labels.csv | 18 |
| Source_Data_ED_Fig8.xlsx | ED Fig. 8b | extfig8b_seaad_tool_transfer.csv | 6 |
| Source_Data_ED_Fig9.xlsx | ED Fig. 9a | extfig9a_cognitive_followup.csv | 4 |
| Source_Data_ED_Fig9.xlsx | ED Fig. 9b | extfig9b_adni_followup_associations.csv | 12 |
| Source_Data_Fig2.xlsx | Fig. 2c | fig2c_doi_internal_pearson.csv | 10 |
| Source_Data_Fig2.xlsx | Fig. 2d | fig2d_doi_internal_error.csv | 20 |
| Source_Data_Fig2.xlsx | Fig. 2g | fig2g_lymphnode_accuracy.csv | 10 |
| Source_Data_Fig2.xlsx | Fig. 2h | fig2h_lymphnode_weightedF1.csv | 10 |
| Source_Data_Fig2.xlsx | Fig. 2i | fig2i_lymphnode_pearson.csv | 10 |
| Source_Data_Fig3.xlsx | Fig. 3c | fig3c_fattyliver_accuracy.csv | 10 |
| Source_Data_Fig3.xlsx | Fig. 3d | fig3d_fattyliver_MCC.csv | 10 |
| Source_Data_Fig3.xlsx | Fig. 3e | fig3e_fattyliver_kappa.csv | 10 |
| Source_Data_Fig3.xlsx | Fig. 3i | fig3i_hypertrophy_accuracy.csv | 10 |
| Source_Data_Fig3.xlsx | Fig. 3j | fig3j_hypertrophy_MCC.csv | 10 |
| Source_Data_Fig3.xlsx | Fig. 3k | fig3k_hypertrophy_kappa.csv | 10 |
| Source_Data_Fig4.xlsx | Fig. 4b | fig4b_glomerulus_counts.csv | 20 |
| Source_Data_Fig5.xlsx | Fig. 5b | fig5bc_cellcounting_overtime.csv | 8 |
| Source_Data_Fig5.xlsx | Fig. 5c | fig5bc_cellcounting_overtime.csv | 7 |
| Source_Data_Fig5.xlsx | Fig. 5f | fig5f_cellcounting_visium.csv | 7 |
| Source_Data_Fig5.xlsx | Fig. 5g | fig5g_cellproportion_visium.csv | 19 |
| Source_Data_Fig5.xlsx | Fig. 5j | fig5jl_cellproportion_overtime.csv | 13 |
| Source_Data_Fig5.xlsx | Fig. 5k | fig5jl_cellproportion_overtime.csv | 31 |
| Source_Data_Fig5.xlsx | Fig. 5l | fig5jl_cellproportion_overtime.csv | 30 |
| Source_Data_Fig6.xlsx | Fig. 6b | fig6b_metastasis_accuracy.csv | 10 |
| Source_Data_Fig6.xlsx | Fig. 6c | fig6c_metastasis_pearson.csv | 10 |
| Source_Data_Fig6.xlsx | Fig. 6d | fig6d_metastasis_radar.csv | 50 |
| Source_Data_Fig6.xlsx | Fig. 6f | fig6f_metastasis_skill_pearson.csv | 5 |
| Source_Data_Fig7.xlsx | Fig. 7d | fig7d_tubule_run1_auc_trajectory.csv | 142 |
| Source_Data_Fig7.xlsx | Fig. 7e | fig7e_tubule_evolution_pearson.csv | 20 |
| Source_Data_Fig7.xlsx | Fig. 7f | fig7f_tubule_evolution_auc.csv | 20 |
| Source_Data_Fig7.xlsx | Fig. 7g | fig7g_tubule_evolution_mae.csv | 20 |
| Source_Data_Fig7.xlsx | Fig. 7h | fig7h_tubule_comparison_auc.csv | 8 |
| Source_Data_Fig7.xlsx | Fig. 7i | fig7i_tubule_comparison_mae.csv | 8 |
| Source_Data_Supp_Fig3.xlsx | Supp. Fig. 3a | suppfig3a_planner_backbone_variation.csv | 20 |
| Source_Data_Supp_Fig3.xlsx | Supp. Fig. 3b | suppfig3b_reflector_backbone_pearson.csv | 10 |
| Source_Data_Supp_Fig3.xlsx | Supp. Fig. 3c | suppfig3c_reflector_vs_noreflector.csv | 4 |
