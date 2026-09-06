set -euo pipefail
python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
old = '''def _top_error_text(best: dict[str, Any], cohort_table: pd.DataFrame) -> str:
    loo = pd.DataFrame(best.get("per_donor_loo", []))
    if loo.empty:
        return "No analyzable donors remained after filtering."
    top = loo.sort_values("abs_error", ascending=False).head(3).merge(
        cohort_table[["donor_id", "cognitive_status", "braak_label", "cerad_label"]],
        on="donor_id",
        how="left",
    )
    parts = []
    for _, row in top.iterrows():
        parts.append(
            f"{row['donor_id']} (abs error {row['abs_error']:.3f}, "
            f"{row.get('cognitive_status', 'NA')}, Braak {row.get('braak_label', 'NA')}, "
            f"CERAD {row.get('cerad_label', 'NA')})"
        )
    return "; ".join(parts)


def _shared_pattern_text(best: dict[str, Any], cohort_table: pd.DataFrame) -> str:
    loo = pd.DataFrame(best.get("per_donor_loo", []))
    if loo.empty:
        return "No clear shared pattern was estimable."
    top = loo.sort_values("abs_error", ascending=False).head(5).merge(
        cohort_table[["donor_id", "cognitive_status", "braak_label", "cerad_label", "overall_ad_neuropath_change"]],
        on="donor_id",
        how="left",
    )
    messages = []
    for col, label in [
        ("cognitive_status", "cognitive status"),
        ("overall_ad_neuropath_change", "AD neuropath change"),
        ("braak_label", "Braak stage"),
        ("cerad_label", "CERAD"),
    ]:
        mode = top[col].mode(dropna=True)
        if not mode.empty:
            frac = float((top[col] == mode.iloc[0]).mean())
            if frac >= 0.6:
                messages.append(f"Most high-error donors share {label}={mode.iloc[0]}")
    if not messages:
        return "The largest errors are donor-specific rather than concentrated in one obvious clinicopathologic stratum."
    return "; ".join(messages) + "."
'''
new = '''def _cohort_meta_subset(cohort_table: pd.DataFrame, preferred: list[str]) -> pd.DataFrame:
    keep = [col for col in preferred if col in cohort_table.columns]
    if "donor_id" not in keep:
        keep = ["donor_id", *keep]
    return cohort_table.loc[:, keep].copy()


def _top_error_text(best: dict[str, Any], cohort_table: pd.DataFrame) -> str:
    loo = pd.DataFrame(best.get("per_donor_loo", []))
    if loo.empty:
        return "No analyzable donors remained after filtering."
    meta = _cohort_meta_subset(cohort_table, ["cognitive_status", "overall_ad_neuropath_change", "braak_numeric", "cerad_ordinal"])
    top = loo.sort_values("abs_error", ascending=False).head(3).merge(meta, on="donor_id", how="left")
    parts = []
    for _, row in top.iterrows():
        details = [f"abs error {row['abs_error']:.3f}"]
        for col, label in [
            ("cognitive_status", "status"),
            ("overall_ad_neuropath_change", "ADNC"),
            ("braak_numeric", "Braak"),
            ("cerad_ordinal", "CERAD"),
        ]:
            if col in row.index and pd.notna(row[col]):
                details.append(f"{label} {row[col]}")
        parts.append(f"{row['donor_id']} ({', '.join(details)})")
    return "; ".join(parts)


def _shared_pattern_text(best: dict[str, Any], cohort_table: pd.DataFrame) -> str:
    loo = pd.DataFrame(best.get("per_donor_loo", []))
    if loo.empty:
        return "No clear shared pattern was estimable."
    meta = _cohort_meta_subset(cohort_table, ["cognitive_status", "overall_ad_neuropath_change", "braak_numeric", "cerad_ordinal"])
    top = loo.sort_values("abs_error", ascending=False).head(5).merge(meta, on="donor_id", how="left")
    messages = []
    for col, label in [
        ("cognitive_status", "cognitive status"),
        ("overall_ad_neuropath_change", "AD neuropath change"),
        ("braak_numeric", "Braak stage"),
        ("cerad_ordinal", "CERAD"),
    ]:
        if col not in top.columns:
            continue
        mode = top[col].mode(dropna=True)
        if not mode.empty:
            frac = float((top[col] == mode.iloc[0]).mean())
            if frac >= 0.6:
                messages.append(f"Most high-error donors share {label}={mode.iloc[0]}")
    if not messages:
        return "The largest errors are donor-specific rather than concentrated in one obvious clinicopathologic stratum."
    return "; ".join(messages) + "."
'''
text = text.replace(old, new)
path.write_text(text)
PY
python /scratch/result.py
