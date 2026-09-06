set -euo pipefail
python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
t=p.read_text()
old="""    dementia_n = int((top_err.get(\"cognitive_status\") == \"Dementia\").sum()) if \"cognitive_status\" in top_err.columns else 0
    high_path_n = (
        int(top_err[\"overall_ad_neuropath_change\"].isin([\"High\", \"Intermediate\"]).sum())
        if \"overall_ad_neuropath_change\" in top_err.columns
        else 0
    )
    severe_col = \"ca1_severe_reactive_cuffed_pyramidal_n_r50um_ra2\"
    severe_med = float(top_err[severe_col].median()) if len(top_err) and severe_col in top_err.columns else float(\"nan\")

    direction = (
"""
new="""    dementia_available = \"cognitive_status\" in top_err.columns and top_err[\"cognitive_status\"].notna().any()
    path_available = \"overall_ad_neuropath_change\" in top_err.columns and top_err[\"overall_ad_neuropath_change\"].notna().any()
    dementia_n = int((top_err[\"cognitive_status\"] == \"Dementia\").sum()) if dementia_available else None
    high_path_n = (
        int(top_err[\"overall_ad_neuropath_change\"].isin([\"High\", \"Intermediate\"]).sum())
        if path_available
        else None
    )
    severe_col = \"ca1_severe_reactive_cuffed_pyramidal_n_r50um_ra2\"
    severe_med = float(top_err[severe_col].median()) if len(top_err) and severe_col in top_err.columns else float(\"nan\")

    if dementia_available and path_available:
        error_shared = (
            f\"Among these five donors, {dementia_n}/5 were Dementia and {high_path_n}/5 had Intermediate-or-High AD neuropath change. \"
        )
    else:
        error_shared = (
            \"The compact cohort table available to this worker does not include cognitive-status or AD-neuropath-change labels, so shared clinical metadata for these outliers could not be checked directly here. \"
        )

    direction = (
"""
t=t.replace(old,new)
old2="""3. Error pattern: Largest absolute LOO errors were in {donor_list}. Among these five donors, {dementia_n}/5 were Dementia and {high_path_n}/5 had Intermediate-or-High AD neuropath change. Their median severe-cuffed CA1 pyramidal count was {_fmt(severe_med)}, suggesting the biomarker still misses some donor-to-donor heterogeneity even when the severe cuff population is abundant.
"""
new2="""3. Error pattern: Largest absolute LOO errors were in {donor_list}. {error_shared}Their median severe-cuffed CA1 pyramidal count was {_fmt(severe_med)}, suggesting the biomarker still misses some donor-to-donor heterogeneity even when the severe cuff population is abundant.
"""
t=t.replace(old2,new2)
p.write_text(t)
PY
python /scratch/result.py >/scratch/run_stdout.txt
sed -n '1,220p' /scratch/report.md
