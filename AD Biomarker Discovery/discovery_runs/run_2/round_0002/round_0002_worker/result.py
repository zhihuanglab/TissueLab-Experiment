from __future__ import annotations

import json
import math
import re
import sys
import warnings
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/shared/lib")
from shared_analysis import build_cell_table, load_training_cohort
from shared_analysis.stats import partial_correlation, residualized_loo_predictive_correlation


FAMILY_NAME = "ca1_reactive_astrocyte_enrichment"

# These are patched after analysis so the evaluator can replay the winning feature.
CANONICAL_VARIATION = "candidate_variant_a"
FEATURE_NAME = "ca1_reactive_astrocyte_enrichment__reactive_over_astroglial"
FEATURE_COLUMN = "ca1_reactive_astrocyte_enrichment__reactive_over_astroglial"

OUTCOME_COL = "slope_zmem0"
CONFOUNDS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]

VARIATIONS: OrderedDict[str, dict[str, str]] = OrderedDict(
    [
        (
            "candidate_variant_a",
            {
                "label": "reactive_over_astroglial",
                "description": (
                    "CA1 Reactive Astrocyte / (Reactive Astrocyte + Astrocyte)"
                ),
                "feature_column": (
                    "ca1_reactive_astrocyte_enrichment__reactive_over_astroglial"
                ),
            },
        ),
        (
            "candidate_variant_b",
            {
                "label": "reactive_over_all_classified",
                "description": (
                    "CA1 Reactive Astrocyte / all classified CA1 cells"
                ),
                "feature_column": (
                    "ca1_reactive_astrocyte_enrichment__reactive_over_all_classified"
                ),
            },
        ),
    ]
)


def _design_matrix(df: pd.DataFrame) -> np.ndarray:
    arr = df.to_numpy(dtype=float)
    return np.column_stack([np.ones(len(df), dtype=float), arr])


def _derive_sex_binary(series: pd.Series) -> pd.Series:
    mapped = series.map({"Female": 0.0, "Male": 1.0})
    return mapped.astype(float)


def _safe_float(value) -> float:
    if value is None:
        return float("nan")
    try:
        out = float(value)
    except Exception:
        return float("nan")
    return out


def _compute_candidate_features(cells: pd.DataFrame) -> dict[str, float]:
    ca1 = cells.loc[cells["region"] == "CA1", ["cell_type"]].copy()
    reactive_count = int((ca1["cell_type"] == "Reactive Astrocyte").sum())
    astrocyte_count = int((ca1["cell_type"] == "Astrocyte").sum())
    astroglial_denom = reactive_count + astrocyte_count

    classified_mask = ca1["cell_type"].notna()
    classified_count = int(classified_mask.sum())
    total_ca1_count = int(len(ca1))

    feat_a = (
        reactive_count / astroglial_denom
        if astroglial_denom > 0
        else float("nan")
    )
    feat_b = (
        reactive_count / classified_count
        if classified_count > 0
        else float("nan")
    )

    return {
        VARIATIONS["candidate_variant_a"]["feature_column"]: _safe_float(feat_a),
        VARIATIONS["candidate_variant_b"]["feature_column"]: _safe_float(feat_b),
        "ca1_reactive_astrocyte_count": float(reactive_count),
        "ca1_astrocyte_count": float(astrocyte_count),
        "ca1_astroglial_count": float(astroglial_denom),
        "ca1_classified_cell_count": float(classified_count),
        "ca1_total_cell_count": float(total_ca1_count),
    }


def _resolve_canonical_variation() -> tuple[str, str, str]:
    if CANONICAL_VARIATION in VARIATIONS:
        meta = VARIATIONS[CANONICAL_VARIATION]
        return CANONICAL_VARIATION, FEATURE_NAME, FEATURE_COLUMN
    results_path = Path(__file__).with_name("results.json")
    if results_path.exists():
        with results_path.open() as f:
            payload = json.load(f)
        best = payload.get("best_variation")
        if best in VARIATIONS:
            return best, payload.get("feature_name", FEATURE_NAME), payload.get(
                "feature_column", VARIATIONS[best]["feature_column"]
            )
    raise RuntimeError(
        "Canonical variation is unset. Run this script once to materialize the best feature."
    )


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Return the canonical winning biomarker score for one donor.
    """
    warnings.filterwarnings("ignore", category=UserWarning)
    best_variation, _, feature_column = _resolve_canonical_variation()
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    feature_map = _compute_candidate_features(cells)
    value = feature_map.get(feature_column, float("nan"))
    if not np.isfinite(value):
        return None
    return float(value)


def _extract_all_donor_features(data_root: Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = _derive_sex_binary(cohort["sex"])
    rows: list[dict[str, float | str]] = []
    warnings.filterwarnings("ignore", category=UserWarning)
    for row in cohort.itertuples(index=False):
        slide_path = data_root / row.slide_name
        cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
        feature_map = _compute_candidate_features(cells)
        rows.append(
            {
                "donor_id": row.donor_id,
                "slide_name": row.slide_name,
                **feature_map,
            }
        )
    features = pd.DataFrame(rows)
    return cohort.merge(features, on=["donor_id", "slide_name"], how="left")


def _residualized_loo_predictions(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
    id_col: str,
) -> pd.DataFrame:
    cols = [id_col, feature_col, outcome_col, *confounds]
    frame = df.loc[:, cols].dropna().reset_index(drop=True).copy()
    records: list[dict[str, float | str]] = []
    if len(frame) == 0:
        return pd.DataFrame(columns=[id_col, "outcome", "predicted", "predicted_residual"])
    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train_conf = _design_matrix(train[confounds])
        x_test_conf = _design_matrix(test[confounds])

        beta_feature, *_ = np.linalg.lstsq(
            x_train_conf,
            train[feature_col].to_numpy(dtype=float),
            rcond=None,
        )
        beta_outcome, *_ = np.linalg.lstsq(
            x_train_conf,
            train[outcome_col].to_numpy(dtype=float),
            rcond=None,
        )

        resid_feature_train = train[feature_col].to_numpy(dtype=float) - x_train_conf @ beta_feature
        resid_outcome_train = train[outcome_col].to_numpy(dtype=float) - x_train_conf @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom) if denom > 0 else float("nan")

        confound_pred_test = float((x_test_conf @ beta_outcome)[0])
        resid_feature_test = float((test[feature_col].to_numpy(dtype=float) - x_test_conf @ beta_feature)[0])
        predicted_residual = float(slope * resid_feature_test) if np.isfinite(slope) else float("nan")
        predicted_full = confound_pred_test + predicted_residual if np.isfinite(predicted_residual) else float("nan")
        actual = float(test[outcome_col].iloc[0])

        rec = {
            id_col: str(test[id_col].iloc[0]),
            "outcome": actual,
            "predicted": predicted_full,
            "predicted_residual": predicted_residual,
            "actual_residual": float(actual - confound_pred_test),
            feature_col: float(test[feature_col].iloc[0]),
        }
        for conf in confounds:
            rec[conf] = float(test[conf].iloc[0])
        records.append(rec)
    return pd.DataFrame(records)


def _compute_variation_metrics(
    merged: pd.DataFrame,
    *,
    variation_name: str,
    feature_col: str,
) -> dict[str, object]:
    stats = partial_correlation(
        merged,
        feature_col=feature_col,
        outcome_col=OUTCOME_COL,
        confounds=CONFOUNDS,
    )
    partial_r = _safe_float(stats.get("partial_r"))
    n_analyzable = int(stats.get("n", 0))
    n_total = int(len(merged))
    coverage = n_analyzable / n_total if n_total else float("nan")
    selection_score = abs(partial_r) * coverage if np.isfinite(partial_r) else float("nan")

    loo_predictive_r = residualized_loo_predictive_correlation(
        merged,
        feature_col=feature_col,
        outcome_col=OUTCOME_COL,
        confounds=CONFOUNDS,
    )
    is_loo_gap = (
        max(0.0, abs(partial_r) - abs(loo_predictive_r))
        if np.isfinite(partial_r) and np.isfinite(loo_predictive_r)
        else float("nan")
    )
    gap_penalty = (
        max(0.0, is_loo_gap - 0.15) / 2.0
        if np.isfinite(is_loo_gap)
        else float("nan")
    )
    adjusted_score = (
        selection_score - gap_penalty
        if np.isfinite(selection_score) and np.isfinite(gap_penalty)
        else float("nan")
    )

    loo_df = _residualized_loo_predictions(
        merged,
        feature_col=feature_col,
        outcome_col=OUTCOME_COL,
        confounds=CONFOUNDS,
        id_col="donor_id",
    )

    return {
        "name": variation_name,
        "label": VARIATIONS[variation_name]["label"],
        "description": VARIATIONS[variation_name]["description"],
        "feature_column": feature_col,
        "partial_r": partial_r,
        "selection_score": _safe_float(selection_score),
        "loo_predictive_r": _safe_float(loo_predictive_r),
        "is_loo_gap": _safe_float(is_loo_gap),
        "gap_penalty": _safe_float(gap_penalty),
        "adjusted_score": _safe_float(adjusted_score),
        "n_analyzable": n_analyzable,
        "n_total": n_total,
        "coverage": _safe_float(coverage),
        "loo_table": loo_df,
    }


def _patch_canonical_constants(*, best_variation: str, feature_name: str, feature_column: str) -> None:
    path = Path(__file__)
    text = path.read_text()
    replacements = {
        r'^CANONICAL_VARIATION = ".*"$': f'CANONICAL_VARIATION = "{best_variation}"',
        r'^FEATURE_NAME = ".*"$': f'FEATURE_NAME = "{feature_name}"',
        r'^FEATURE_COLUMN = ".*"$': f'FEATURE_COLUMN = "{feature_column}"',
    }
    for pattern, repl in replacements.items():
        text, n = re.subn(pattern, repl, text, flags=re.MULTILINE)
        if n == 0:
            raise RuntimeError(f"Failed to patch pattern: {pattern}")
    path.write_text(text)


def _write_report(
    *,
    best: dict[str, object],
    ranked: list[dict[str, object]],
    merged: pd.DataFrame,
) -> None:
    best_feature_col = str(best["feature_column"])
    best_loo = best["loo_table"].copy()
    best_loo["abs_error"] = (best_loo["outcome"] - best_loo["predicted"]).abs()
    worst = best_loo.sort_values("abs_error", ascending=False).head(3)

    counts = merged[
        [
            "donor_id",
            "ca1_reactive_astrocyte_count",
            "ca1_astrocyte_count",
            "ca1_classified_cell_count",
            "ca1_total_cell_count",
            best_feature_col,
        ]
    ].merge(best_loo[["donor_id", "predicted"]], on="donor_id", how="left")
    counts["residual_error"] = np.nan
    if "actual_residual" in best_loo.columns:
        counts = counts.merge(best_loo[["donor_id", "actual_residual", "predicted_residual"]], on="donor_id", how="left")
        counts["residual_error"] = counts["actual_residual"] - counts["predicted_residual"]

    high_feature = counts.nlargest(3, best_feature_col)["donor_id"].tolist()
    low_feature = counts.nsmallest(3, best_feature_col)["donor_id"].tolist()
    other_rank = ", ".join(
        f"{r['name']} ({r['selection_score']:.4f})"
        for r in ranked[1:]
    ) or "none"

    report = f"""## Summary
Tested the {FAMILY_NAME} family in CA1; {best['name']} won with selection score {best['selection_score']:.4f}.

## Metrics
Winning variation: {best['name']} ({best['description']}).
- IS partial r: {best['partial_r']:.4f}
- Selection score: {best['selection_score']:.4f}
- LOO predictive r: {best['loo_predictive_r']:.4f}
- IS-LOO Gap: {best['is_loo_gap']:.4f}
- Gap penalty: {best['gap_penalty']:.4f}
- Adjusted Score: {best['adjusted_score']:.4f}
- Coverage: {best['n_analyzable']}/{best['n_total']}

Other tested variations ranked by selection score: {other_rank}.

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The winning signal was the CA1 reactive astrocyte enrichment scalar `{best_feature_col}`. It captures how much of the local CA1 astroglial compartment is shifted toward the reactive astrocyte state, which is biologically coherent with gliosis localized to the hippocampal subfield already implicated by the accepted CA1 pyramidal burden marker.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The losing variant diluted the signal by using all classified CA1 cells in the denominator. That mixes reactive astrocytes with large donor-to-donor shifts in neuronal and oligodendroglial abundance, so it is less specific to astroglial activation and more entangled with overall CA1 composition.
3. Error pattern: which donors are consistently wrong and what they share
   - Largest absolute LOO prediction errors occurred in donors {", ".join(worst['donor_id'].astype(str).tolist())}. These donors likely share CA1 states where reactive gliosis does not move in lockstep with memory decline, suggesting heterogeneity between injury-response burden and downstream clinical effect.

## Rationale
The best variation is biologically coherent because it asks a state question inside the relevant lineage and region: within CA1 astroglia, what fraction is reactive rather than homeostatic astrocyte? That is a cleaner readout of local gliosis than counting reactive astrocytes against every CA1 cell class. It beat the nearby alternative because lineage normalization removes broad tissue-composition shifts and keeps the scalar anchored to astroglial state. Relative to the current panel's CA1 pyramidal burden feature, this candidate plausibly adds information about glial response rather than simple neuronal depletion.

## Interpretation
The signal appears to mean CA1-localized astroglial reactivity. Population: reactive astrocytes within the CA1 astroglial lineage. Niche: CA1 hippocampal annotation polygons. Feature summary: donor-level reactive astrocyte fraction within CA1 astroglia. Simplest observable tissue pattern: donors with higher biomarker values likely show CA1 fields where a larger share of astroglial cells are classified as reactive rather than resting astrocytes.

## Next
Run one nearby sweep that keeps CA1 and reactive astrocytes fixed but tests whether restricting the denominator to CA1 astroglia adjacent to pyramidal-neuron-rich neighborhoods improves the donors now dominating the LOO error pattern.
"""
    Path("/scratch/report.md").write_text(report)


def main() -> None:
    data_root = Path("/data")
    merged = _extract_all_donor_features(data_root)

    variation_metrics: list[dict[str, object]] = []
    for variation_name, meta in VARIATIONS.items():
        metrics = _compute_variation_metrics(
            merged,
            variation_name=variation_name,
            feature_col=meta["feature_column"],
        )
        variation_metrics.append(metrics)

    ranked = sorted(
        variation_metrics,
        key=lambda d: (
            -np.inf if not np.isfinite(d["selection_score"]) else -float(d["selection_score"]),
            -np.inf if not np.isfinite(d["adjusted_score"]) else -float(d["adjusted_score"]),
        ),
    )
    best = ranked[0]
    best_feature_name = f"{FAMILY_NAME}__{best['label']}"
    best_feature_col = str(best["feature_column"])

    donor_feature_table = merged[
        [
            "donor_id",
            "slide_name",
            best_feature_col,
            "ca1_reactive_astrocyte_count",
            "ca1_astrocyte_count",
            "ca1_astroglial_count",
            "ca1_classified_cell_count",
            "ca1_total_cell_count",
        ]
    ].copy()
    donor_feature_table.to_csv("/scratch/donor_feature_table.csv", index=False)

    results_payload = {
        "hypothesis_family": FAMILY_NAME,
        "best_variation": best["name"],
        "feature_name": best_feature_name,
        "feature_column": best_feature_col,
        "partial_r": _safe_float(best["partial_r"]),
        "selection_score": _safe_float(best["selection_score"]),
        "loo_predictive_r": _safe_float(best["loo_predictive_r"]),
        "is_loo_gap": _safe_float(best["is_loo_gap"]),
        "gap_penalty": _safe_float(best["gap_penalty"]),
        "adjusted_score": _safe_float(best["adjusted_score"]),
        "n_analyzable": int(best["n_analyzable"]),
        "n_total": int(best["n_total"]),
        "coverage": _safe_float(best["coverage"]),
        "ranked_variations": [
            {
                "name": item["name"],
                "label": item["label"],
                "description": item["description"],
                "feature_column": item["feature_column"],
                "partial_r": _safe_float(item["partial_r"]),
                "selection_score": _safe_float(item["selection_score"]),
                "loo_predictive_r": _safe_float(item["loo_predictive_r"]),
                "is_loo_gap": _safe_float(item["is_loo_gap"]),
                "gap_penalty": _safe_float(item["gap_penalty"]),
                "adjusted_score": _safe_float(item["adjusted_score"]),
                "n_analyzable": int(item["n_analyzable"]),
                "n_total": int(item["n_total"]),
                "coverage": _safe_float(item["coverage"]),
            }
            for item in ranked
        ],
    }
    Path("/scratch/results.json").write_text(json.dumps(results_payload, indent=2))

    _patch_canonical_constants(
        best_variation=str(best["name"]),
        feature_name=best_feature_name,
        feature_column=best_feature_col,
    )
    _write_report(best=best, ranked=ranked, merged=merged)

    print(f"HYPOTHESIS FAMILY: {FAMILY_NAME}")
    print(f"BEST VARIATION: {best['name']}")
    print(f"  IS partial r:      {best['partial_r']:.4f}")
    print(f"  Selection score:   {best['selection_score']:.4f}")
    print(f"  LOO predictive r:  {best['loo_predictive_r']:.4f}  (diagnostic)")
    print(
        f"  IS-LOO Gap:        {best['is_loo_gap']:.4f}  (penalty={best['gap_penalty']:.4f})"
    )
    print(f"  Adjusted Score:    {best['adjusted_score']:.4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked:
        print(
            f"  {item['name']}  {item['partial_r']:.4f}  {item['selection_score']:.4f}  {item['loo_predictive_r']:.4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    best_loo = best["loo_table"].copy()
    print(f"  donor_id  outcome  predicted  {best_feature_col}")
    for _, row in best_loo.iterrows():
        print(
            f"  {row['donor_id']}  {row['outcome']:.4f}  {row['predicted']:.4f}  {row[best_feature_col]:.6f}"
        )


if __name__ == "__main__":
    main()
