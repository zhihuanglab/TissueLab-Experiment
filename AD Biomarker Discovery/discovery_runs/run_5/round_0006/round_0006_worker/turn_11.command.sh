python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
old = '''def make_error_pattern_text(best_loo: pd.DataFrame) -> str:
    if best_loo.empty:
        return "No analyzable donors for leave-one-out diagnostics."
    err = best_loo.sort_values("absolute_error", ascending=False).head(3)
    donors = []
    for _, row in err.iterrows():
        donors.append(
            f'{row["donor_id"]} (abs err {row["absolute_error"]:.3f}, '
            f'braak {int(row["braak_numeric"])}, cerad {int(row["cerad_ordinal"])}, '
            f'{row["cognitive_status"]})'
        )
    statuses = err["cognitive_status"].value_counts().to_dict()
    pathology = err["overall_ad_neuropath_change"].value_counts().to_dict()
    return (
        "Largest LOO misses were "
        + "; ".join(donors)
        + f". Shared pattern: cognitive-status counts {statuses} and AD-neuropath-change counts {pathology}."
    )
'''
new = '''def make_error_pattern_text(best_loo: pd.DataFrame) -> str:
    if best_loo.empty:
        return "No analyzable donors for leave-one-out diagnostics."
    err = best_loo.sort_values("absolute_error", ascending=False).head(3)
    donors = []
    for _, row in err.iterrows():
        braak = row.get("braak_numeric")
        cerad = row.get("cerad_ordinal")
        cog = row.get("cognitive_status", "unknown")
        neuropath = row.get("overall_ad_neuropath_change", "unknown")
        donors.append(
            f'{row["donor_id"]} (abs err {row["absolute_error"]:.3f}, '
            f'braak {int(braak) if pd.notna(braak) else "NA"}, '
            f'cerad {int(cerad) if pd.notna(cerad) else "NA"}, '
            f'{cog}, {neuropath})'
        )
    status_col = "cognitive_status" if "cognitive_status" in err.columns else None
    path_col = "overall_ad_neuropath_change" if "overall_ad_neuropath_change" in err.columns else None
    statuses = err[status_col].value_counts().to_dict() if status_col else {}
    pathology = err[path_col].value_counts().to_dict() if path_col else {}
    return (
        "Largest LOO misses were "
        + "; ".join(donors)
        + f". Shared pattern: cognitive-status counts {statuses} and AD-neuropath-change counts {pathology}."
    )
'''
text = text.replace(old, new)
path.write_text(text)
PY
python /scratch/result.py
echo '--- report.md ---'
sed -n '1,240p' /scratch/report.md
