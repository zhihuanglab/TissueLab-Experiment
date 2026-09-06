set -euo pipefail
python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
text=p.read_text()
text=text.replace(
"""    meta_cols = [
        \"donor_id\",
        \"cognitive_status\",
        \"overall_ad_neuropath_change\",
        \"braak_numeric\",
        \"cerad_ordinal\",
        \"ca1_severe_reactive_cuffed_pyramidal_n_r50um_ra2\",
        \"ca1_conditional_immune_cuff_positive_n_r50um_ra2_ly1\",
        \"ca1_conditional_immune_cuff_positive_n_r50um_ra2_ly2\",
    ]
    meta = donor_df.reindex(columns=meta_cols)
    merged = predictions_table.merge(meta, on=\"donor_id\", how=\"left\")
""",
"""    meta_cols = [\"donor_id\", \"cognitive_status\", \"overall_ad_neuropath_change\", \"braak_numeric\", \"cerad_ordinal\"]
    meta = donor_df.reindex(columns=meta_cols)
    merged = predictions_table.merge(meta, on=\"donor_id\", how=\"left\")
""")
text=text.replace(
"""    dementia_n = int((top_err[\"cognitive_status\"] == \"Dementia\").sum())
    high_path_n = int(top_err[\"overall_ad_neuropath_change\"].isin([\"High\", \"Intermediate\"]).sum())
    severe_med = float(top_err[\"ca1_severe_reactive_cuffed_pyramidal_n_r50um_ra2\"].median()) if len(top_err) else float(\"nan\")
""",
"""    dementia_n = int((top_err.get(\"cognitive_status\") == \"Dementia\").sum()) if \"cognitive_status\" in top_err.columns else 0
    high_path_n = (
        int(top_err[\"overall_ad_neuropath_change\"].isin([\"High\", \"Intermediate\"]).sum())
        if \"overall_ad_neuropath_change\" in top_err.columns
        else 0
    )
    severe_col = \"ca1_severe_reactive_cuffed_pyramidal_n_r50um_ra2\"
    severe_med = float(top_err[severe_col].median()) if len(top_err) and severe_col in top_err.columns else float(\"nan\")
""")
p.write_text(text)
PY
python /scratch/result.py
ls -l /scratch/result.py /scratch/results.json /scratch/report.md /scratch/donor_feature_table.csv /scratch/best_variation.json
