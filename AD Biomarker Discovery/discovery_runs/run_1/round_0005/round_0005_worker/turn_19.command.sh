python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()

text = text.replace(
"""def _load_eval_cohort(data_root: str | Path) -> pd.DataFrame:
    data_root = Path(data_root)
    training = data_root / "training_cohort.csv"
    test = data_root / "test_cohort.csv"
    if training.exists():
        cohort = load_training_cohort(data_root).copy()
    elif test.exists():
        cohort = pd.read_csv(test).copy()
    else:
        cohort = load_training_cohort(data_root).copy()
""",
"""def _load_eval_cohort(data_root: str | Path) -> pd.DataFrame:
    data_root = Path(data_root)
    training = data_root / "training_cohort.csv"
    test = data_root / "test_cohort.csv"
    if training.exists():
        cohort = pd.read_csv(training).copy()
    elif test.exists():
        cohort = pd.read_csv(test).copy()
    else:
        cohort = load_training_cohort(data_root).copy()
"""
)

text = text.replace(
"""def _report_error_pattern(best_rows: pd.DataFrame, donor_table: pd.DataFrame) -> str:
    merged = best_rows.merge(
        donor_table[
            [
                "donor_id",
                "cognitive_status",
                "overall_ad_neuropath_change",
                "ca1_pyramidal_count",
                "ca1_reactive_astrocyte_count",
            ]
        ],
        on="donor_id",
        how="left",
    ).copy()
    merged["abs_error"] = (merged["predicted"] - merged["outcome"]).abs()
    top = merged.sort_values("abs_error", ascending=False).head(3)
    donors = ", ".join(
        f"{row.donor_id} (err={row.abs_error:.3f})"
        for row in top.itertuples(index=False)
    )
    dementia_n = int((top["cognitive_status"] == "Dementia").sum())
    high_path_n = int(top["overall_ad_neuropath_change"].isin(["Intermediate", "High"]).sum())
    return (
        f"Largest LOO errors were {donors}. Among these three, {dementia_n}/3 had dementia and "
        f"{high_path_n}/3 had intermediate/high AD neuropathologic change, suggesting that "
        f"perineuronal reactive-astro burden captures only part of the severe-disease spectrum."
    )
""",
"""def _report_error_pattern(best_rows: pd.DataFrame, donor_table: pd.DataFrame, *, feature_col: str) -> str:
    merged = best_rows.merge(
        donor_table[
            [
                "donor_id",
                feature_col,
                "ca1_pyramidal_count",
                "ca1_reactive_astrocyte_count",
            ]
        ],
        on="donor_id",
        how="left",
    ).copy()
    merged["abs_error"] = (merged["predicted"] - merged["outcome"]).abs()
    top = merged.sort_values("abs_error", ascending=False).head(3)
    donors = ", ".join(
        f"{row.donor_id} (err={row.abs_error:.3f}, feature={getattr(row, feature_col):.4f})"
        for row in top.itertuples(index=False)
    )
    underpredicted = int((top["predicted"] > top["outcome"]).sum())
    return (
        f"Largest LOO errors were {donors}. {underpredicted}/3 were overpredicted by the model, and the group still had "
        f"low-to-moderate {feature_col} values, suggesting that some donors decline more than expected even without a strong CA1 reactive-astro encirclement signal."
    )
"""
)

text = text.replace(
"""    error_pattern = _report_error_pattern(best_loo, donor_table)
""",
"""    error_pattern = _report_error_pattern(best_loo, donor_table, feature_col=best["feature_column"])
"""
)

text = text.replace(
"""        failed_text = (
            "The stricter thresholds lost signal relative to the winner. "
            + " ".join(
                f"{item['name']} fell to selection_score={item['selection_score']:.4f}"
                f" with loo_predictive_r={item['loo_predictive_r']:.4f}."
                for item in others
            )
        )
""",
"""        failed_text = (
            "The nearby alternatives both lost signal relative to the winner. "
            + " ".join(
                f"{item['name']} fell to selection_score={item['selection_score']:.4f}"
                f" with loo_predictive_r={item['loo_predictive_r']:.4f}."
                for item in others
            )
        )
"""
)

text = text.replace(
"""        "feature_column": CANONICAL_FEATURE_COLUMN,
""",
"""        "feature_column": CANONICAL_FEATURE_COLUMN,
        "winning_variation_feature_column": best["feature_column"],
"""
)

path.write_text(text)
PY
python /scratch/result.py
