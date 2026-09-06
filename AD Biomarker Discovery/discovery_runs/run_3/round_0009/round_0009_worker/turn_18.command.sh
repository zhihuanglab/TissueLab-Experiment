set -euo pipefail
python - <<'PY'
from pathlib import Path
p = Path('/scratch/result.py')
text = p.read_text()
old = """    merged = predictions_table.merge(
        donor_df[
            [
                \"donor_id\",
                \"cognitive_status\",
                \"overall_ad_neuropath_change\",
                \"braak_numeric\",
                \"cerad_ordinal\",
                \"ca1_severe_reactive_cuffed_pyramidal_n_r50um_ra2\",
                \"ca1_conditional_immune_cuff_positive_n_r50um_ra2_ly1\",
                \"ca1_conditional_immune_cuff_positive_n_r50um_ra2_ly2\",
            ]
        ],
        on=\"donor_id\",
        how=\"left\",
    )
"""
new = """    meta_cols = [
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
"""
p.write_text(text.replace(old, new))
PY
python /scratch/result.py
ls -l /scratch/result.py /scratch/results.json /scratch/report.md /scratch/donor_feature_table.csv /scratch/best_variation.json
