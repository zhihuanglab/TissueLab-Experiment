set -euo pipefail
python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
txt=p.read_text()
old = """    merged = predictions_table.merge(
        donor_df[
            [
                \"donor_id\",
                \"cognitive_status\",
                \"overall_ad_neuropath_change\",
                \"braak_numeric\",
                \"cerad_ordinal\",
                \"ca1_severe_reactive_cuffed_pyramidal_n_r50um_ra2\",
                \"ca1_conditional_immune_cuff_positive_n_r50um_ra2_ly2\",
                \"ca1_conditional_immune_cuff_positive_n_r50um_ra2_ly3\",
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
        \"ca1_conditional_immune_cuff_positive_n_r50um_ra2_ly2\",
        \"ca1_conditional_immune_cuff_positive_n_r50um_ra2_ly3\",
    ]
    merged = predictions_table.merge(
        donor_df.reindex(columns=meta_cols),
        on=\"donor_id\",
        how=\"left\",
    )
"""
txt=txt.replace(old,new)
old2 = """    dementia_n = int((top_err[\"cognitive_status\"] == \"Dementia\").sum())
    high_path_n = int(top_err[\"overall_ad_neuropath_change\"].isin([\"High\", \"Intermediate\"]).sum())
"""
new2 = """    dementia_available = \"cognitive_status\" in top_err.columns and top_err[\"cognitive_status\"].notna().any()
    path_available = \"overall_ad_neuropath_change\" in top_err.columns and top_err[\"overall_ad_neuropath_change\"].notna().any()
    dementia_n = int((top_err[\"cognitive_status\"] == \"Dementia\").sum()) if dementia_available else None
    high_path_n = int(top_err[\"overall_ad_neuropath_change\"].isin([\"High\", \"Intermediate\"]).sum()) if path_available else None
"""
txt=txt.replace(old2,new2)
old3 = """3. Error pattern: Largest absolute LOO errors were in {donor_list}. Among these five donors, {dementia_n}/5 were Dementia and {high_path_n}/5 had Intermediate-or-High AD neuropath change. Their median severe-cuffed CA1 pyramidal count was {_fmt(severe_med)}, with median ly>=2-positive count {_fmt(ly2_med)} and ly>=3-positive count {_fmt(ly3_med)}, suggesting residual error is partly driven by donor-to-donor differences in how often the severe cuff progresses into a very dense immune pocket.
"""
new3 = """3. Error pattern: Largest absolute LOO errors were in {donor_list}. {('Among these five donors, ' + str(dementia_n) + '/5 were Dementia and ' + str(high_path_n) + '/5 had Intermediate-or-High AD neuropath change. ') if dementia_available and path_available else 'Clinical-status and AD-neuropath-change labels were not available in the local cohort table used by this worker, so those shared metadata could not be summarized directly here. '}Their median severe-cuffed CA1 pyramidal count was {_fmt(severe_med)}, with median ly>=2-positive count {_fmt(ly2_med)} and ly>=3-positive count {_fmt(ly3_med)}, suggesting residual error is partly driven by donor-to-donor differences in how often the severe cuff progresses into a very dense immune pocket.
"""
txt=txt.replace(old3,new3)
p.write_text(txt)
print('patched')
PY

python /scratch/result.py
printf '\n--- results.json head ---\n'
sed -n '1,160p' /scratch/results.json
printf '\n--- report.md ---\n'
sed -n '1,220p' /scratch/report.md
