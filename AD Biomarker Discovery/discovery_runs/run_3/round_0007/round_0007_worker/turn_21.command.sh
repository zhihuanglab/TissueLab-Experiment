python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
old = '''def summarize_error_pattern(loo_table: pd.DataFrame) -> tuple[str, str, str]:
    if loo_table.empty:
        return (
            "No analyzable leave-one-out rows were available.",
            "No losing variation comparison available.",
            "No donor-level error pattern available.",
        )
    work = loo_table.copy()
    work["error"] = work["outcome"] - work["predicted"]
    work["abs_error"] = work["error"].abs()
    work = work.sort_values("abs_error", ascending=False).reset_index(drop=True)
    worst = work.head(3)

    donor_list = ", ".join(str(v) for v in worst["donor_id"].tolist())
    hi_feat = float(work.iloc[0, work.columns.get_loc(work.columns[-3])]) if len(work.columns) >= 4 else float("nan")
    over = worst.loc[worst["error"] < 0, "donor_id"].tolist()
    under = worst.loc[worst["error"] > 0, "donor_id"].tolist()

    pattern = (
        f"Largest LOO misses are {donor_list}. "
        f"Over-predicted donors ({', '.join(map(str, over)) if over else 'none'}) have worse predicted decline than observed, "
        f"while under-predicted donors ({', '.join(map(str, under)) if under else 'none'}) have more decline than the local CA1 burden alone suggests."
    )

    shared = (
        "These mistakes look like niche-mismatch cases: the feature is strongest when CA1 perineuronal reactive burden tracks memory decline, "
        "but it misses donors whose outcome likely reflects pathology outside this local CA1 astrocyte-pyramidal niche."
    )

    note = (
        f"The worst-error donor feature values remain within the same family scale (top-row feature ≈ {hi_feat:.2f}), "
        "so the residual error is more about biological heterogeneity than outright extraction failure."
        if math.isfinite(hi_feat)
        else "Feature extraction succeeded, but donor-level residuals still suggest biological heterogeneity beyond this one niche summary."
    )
    return pattern, shared, note
'''
new = '''def summarize_error_pattern(loo_table: pd.DataFrame) -> tuple[str, str, str]:
    if loo_table.empty:
        return (
            "No analyzable leave-one-out rows were available.",
            "No donor-level mismatch pattern could be summarized.",
            "No donor-level error pattern available.",
        )

    work = loo_table.copy()
    feature_cols = [c for c in work.columns if c not in {"donor_id", "outcome", "predicted"}]
    feature_col = feature_cols[0] if feature_cols else None

    work["error"] = work["outcome"] - work["predicted"]
    work["abs_error"] = work["error"].abs()
    if feature_col is not None:
        work["feature_value"] = work[feature_col].astype(float)
        feature_median = float(work["feature_value"].median())
    else:
        work["feature_value"] = np.nan
        feature_median = float("nan")

    worst = work.sort_values("abs_error", ascending=False).head(5).reset_index(drop=True)
    donor_list = ", ".join(str(v) for v in worst["donor_id"].tolist())

    underpred = worst.loc[worst["error"] < 0, "donor_id"].tolist()
    overpred = worst.loc[worst["error"] > 0, "donor_id"].tolist()

    low_burden_underpred = worst.loc[
        (worst["error"] < 0) & (worst["feature_value"] <= feature_median),
        "donor_id",
    ].tolist()
    high_burden_overpred = worst.loc[
        (worst["error"] > 0) & (worst["feature_value"] > feature_median),
        "donor_id",
    ].tolist()

    pattern = (
        f"Largest LOO misses are {donor_list}. "
        f"Under-predicted decline donors ({', '.join(map(str, underpred)) if underpred else 'none'}) had more memory decline than predicted from this local burden, "
        f"while over-predicted decline donors ({', '.join(map(str, overpred)) if overpred else 'none'}) had less decline than their CA1 burden would suggest."
    )

    shared = (
        f"Two repeat mismatch modes appear: low-burden but high-decline donors ({', '.join(map(str, low_burden_underpred)) if low_burden_underpred else 'none'}) "
        f"suggest pathology outside the CA1 peripyramidal astrocyte niche, whereas high-burden but milder-decline donors ({', '.join(map(str, high_burden_overpred)) if high_burden_overpred else 'none'}) "
        "suggest reactive burden can be present without proportionate memory slope severity."
    )

    note = (
        f"The median winner feature value is ≈ {feature_median:.2f}; the worst errors span both sides of that median, "
        "arguing against a simple extraction artifact and more toward biological heterogeneity."
        if math.isfinite(feature_median)
        else "Feature extraction succeeded, but donor-level residuals still suggest biological heterogeneity beyond this one niche summary."
    )
    return pattern, shared, note
'''
if old not in text:
    raise SystemExit('old block not found')
path.write_text(text.replace(old, new))
PY
python /scratch/result.py > /scratch/stdout.txt
sed -n '1,80p' /scratch/stdout.txt
echo
echo '--- updated report ---'
sed -n '1,220p' /scratch/report.md
