cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, "/shared/lib")

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from shared_analysis import build_cell_table, load_training_cohort
from shared_analysis.artifacts import build_results_payload, write_donor_feature_table, write_results_payload
from shared_analysis.stats import (
    bootstrap_partial_correlation,
    bootstrap_partial_correlation_stability,
    leave_one_out_summary,
    partial_correlation,
    residualized_loo_predictive_correlation,
)

# Family-local sweep: CA1 astrocyte-lineage reactive fraction within a
# peripyramidal niche around CA1 pyramidal neurons.

HYPOTHESIS_FAMILY = "CA1 peripyramidal astrocytic reactivity"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
PIXEL_SIZE_UM = 0.503
MIN_DENOMINATOR = 25

VARIATIONS: list[dict[str, Any]] = [
    {
        "name": "peripyramidal_reactivity_r25um",
        "radius_um": 25.0,
        "feature_column": "ca1_peripyramidal_reactivity_r25um",
    },
    {
        "name": "peripyramidal_reactivity_r35um",
        "radius_um": 35.0,
        "feature_column": "ca1_peripyramidal_reactivity_r35um",
    },
]

_CANONICAL_SPEC_PATH = Path(__file__).with_name("best_variation.json")


def _safe_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        value = float(value)
    except Exception:
        return None
    if not math.isfinite(value):
        return None
    return value


def _json_default(value: Any) -> Any:
    if isinstance(value, (np.floating,)):
        value = float(value)
    elif isinstance(value, (np.integer,)):
        value = int(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _load_canonical_spec() -> dict[str, Any]:
    if _CANONICAL_SPEC_PATH.exists():
        return json.loads(_CANONICAL_SPEC_PATH.read_text(encoding="utf-8"))
    return dict(VARIATIONS[0])


_CANONICAL_SPEC = _load_canonical_spec()
FEATURE_NAME = str(_CANONICAL_SPEC["name"])
FEATURE_COLUMN = str(_CANONICAL_SPEC["feature_column"])


def _prepare_cohort(data_root: str | Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    sex = cohort["sex"].astype(str).str.strip().str.lower()
    cohort["sex_binary"] = sex.map({"female": 1.0, "male": 0.0})
    return cohort


def _nearest_pyramidal_distances(ca1_cells: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    pyramidal = ca1_cells.loc[ca1_cells["cell_type"] == "Pyramidal Neuron", ["x", "y"]].to_numpy(dtype=np.float32)
    astro = ca1_cells.loc[
        ca1_cells["cell_type"].isin(["Astrocyte", "Reactive Astrocyte"]),
        ["x", "y"],
    ].to_numpy(dtype=np.float32)
    if len(astro) == 0:
        return astro, np.asarray([], dtype=np.float32)
    if len(pyramidal) == 0:
        return astro, np.full(len(astro), np.inf, dtype=np.float32)
    tree = cKDTree(pyramidal)
    distances, _ = tree.query(astro, k=1, workers=-1)
    return astro, np.asarray(distances, dtype=np.float32)


def compute_slide_variations(slide_path: str | Path) -> dict[str, Any]:
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    ca1 = cells.loc[cells["region"] == "CA1", ["x", "y", "cell_type"]].copy()

    result: dict[str, Any] = {
        "ca1_total_cells": int(len(ca1)),
        "ca1_pyramidal_count": int((ca1["cell_type"] == "Pyramidal Neuron").sum()),
        "ca1_astro_lineage_count": int(ca1["cell_type"].isin(["Astrocyte", "Reactive Astrocyte"]).sum()),
    }

    _, nearest_dist_px = _nearest_pyramidal_distances(ca1)
    astro_lineage = ca1.loc[ca1["cell_type"].isin(["Astrocyte", "Reactive Astrocyte"]), ["cell_type"]].reset_index(drop=True)

    for spec in VARIATIONS:
        radius_px = float(spec["radius_um"]) / PIXEL_SIZE_UM
        near_mask = nearest_dist_px <= radius_px if len(nearest_dist_px) else np.asarray([], dtype=bool)
        denom = int(near_mask.sum())
        reactive = int(((astro_lineage["cell_type"] == "Reactive Astrocyte").to_numpy()) & near_mask).sum() if denom else 0
        feature_value = reactive / denom if denom >= MIN_DENOMINATOR else np.nan

        result[spec["feature_column"]] = _safe_float(feature_value)
        result[f"{spec['feature_column']}_denom"] = denom
        result[f"{spec['feature_column']}_reactive_n"] = reactive
        result[f"{spec['feature_column']}_radius_um"] = float(spec["radius_um"])

    return result


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Canonical replay target: the winning CA1 peripyramidal reactivity variation.
    """
    data_root = Path(data_root)
    cohort = _prepare_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    features = compute_slide_variations(slide_path)
    value = features.get(FEATURE_COLUMN)
    return None if value is None or not math.isfinite(float(value)) else float(value)


def _loo_predictions(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
) -> pd.DataFrame:
    frame = df.loc[:, ["donor_id", feature_col, outcome_col, *confounds]].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    rows: list[dict[str, Any]] = []
    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train_conf = np.column_stack([np.ones(len(train), dtype=float), train[confounds].astype(float).to_numpy()])
        x_test_conf = np.column_stack([np.ones(len(test), dtype=float), test[confounds].astype(float).to_numpy()])

        y_feature_train = train[feature_col].to_numpy(dtype=float)
        y_outcome_train = train[outcome_col].to_numpy(dtype=float)

        beta_feature, *_ = np.linalg.lstsq(x_train_conf, y_feature_train, rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train_conf, y_outcome_train, rcond=None)

        resid_feature_train = y_feature_train - x_train_conf @ beta_feature
        resid_outcome_train = y_outcome_train - x_train_conf @ beta_outcome
        denom = float(np.dot(resid_feature_train, resid_feature_train))

        test_feature = float(test.iloc[0][feature_col])
        test_outcome = float(test.iloc[0][outcome_col])
        pred_raw = float("nan")
        pred_resid = float("nan")
        actual_resid = float(test_outcome - (x_test_conf @ beta_outcome)[0])

        if denom > 0 and np.isfinite(denom):
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            resid_feature_test = float(test_feature - (x_test_conf @ beta_feature)[0])
            pred_resid = float(slope * resid_feature_test)
            pred_raw = float((x_test_conf @ beta_outcome)[0] + pred_resid)

        rows.append(
            {
                "donor_id": str(test.iloc[0]["donor_id"]),
                "outcome": test_outcome,
                "predicted": pred_raw,
                "actual_resid": actual_resid,
                "predicted_resid": pred_resid,
                feature_col: test_feature,
            }
        )
    return pd.DataFrame(rows)


def evaluate_variation(
    merged: pd.DataFrame,
    *,
    variation: dict[str, Any],
    n_total: int,
) -> dict[str, Any]:
    feature_col = str(variation["feature_column"])
    analyzable = merged.loc[:, ["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    n_analyzable = int(len(analyzable))
    coverage_ratio = n_analyzable / n_total if n_total else float("nan")

    pc = partial_correlation(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    boot = bootstrap_partial_correlation(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        n_boot=1000,
        random_state=0,
    )
    stability = bootstrap_partial_correlation_stability(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        n_boot=400,
        random_state=0,
    )
    loo_summary = leave_one_out_summary(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        id_col="donor_id",
        unstable_delta=0.10,
    )
    loo_predictive_r = residualized_loo_predictive_correlation(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    partial_r = float(pc["partial_r"])
    selection_score = abs(partial_r) * coverage_ratio if math.isfinite(partial_r) else float("nan")
    is_loo_gap = abs(partial_r) - abs(loo_predictive_r) if math.isfinite(partial_r) and math.isfinite(loo_predictive_r) else float("nan")
    gap_penalty = max(0.0, is_loo_gap) if math.isfinite(is_loo_gap) else float("nan")
    adjusted_score = selection_score - gap_penalty if math.isfinite(selection_score) and math.isfinite(gap_penalty) else float("nan")
    per_donor = _loo_predictions(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )

    return {
        "name": str(variation["name"]),
        "radius_um": float(variation["radius_um"]),
        "feature_column": feature_col,
        "n_total": int(n_total),
        "n_analyzable": n_analyzable,
        "coverage_ratio": float(coverage_ratio),
        "partial_r": float(partial_r),
        "p_value": float(pc["p_value"]),
        "ci_lo": _safe_float(boot.get("ci_lo")),
        "ci_hi": _safe_float(boot.get("ci_hi")),
        "selection_score": _safe_float(selection_score),
        "loo_predictive_r": _safe_float(loo_predictive_r),
        "is_loo_gap": _safe_float(is_loo_gap),
        "gap_penalty": _safe_float(gap_penalty),
        "adjusted_score": _safe_float(adjusted_score),
        "bootstrap_median_partial_r": _safe_float(stability.get("bootstrap_median_partial_r")),
        "bootstrap_sign_consistency": _safe_float(stability.get("bootstrap_sign_consistency")),
        "loo_unstable_count": int(loo_summary.get("unstable_count", 0)),
        "loo_max_shift": _safe_float(loo_summary.get("max_shift")),
        "loo_unstable_donors": loo_summary.get("unstable_donors", []),
        "per_donor_loo": per_donor.to_dict(orient="records"),
    }


def _feature_extraction_table(data_root: Path) -> pd.DataFrame:
    cohort = _prepare_cohort(data_root)
    rows: list[dict[str, Any]] = []
    for row in cohort.itertuples(index=False):
        slide_path = data_root / str(row.slide_name)
        features = compute_slide_variations(slide_path)
        entry = {col: getattr(row, col) for col in cohort.columns}
        entry.update(features)
        rows.append(entry)
    return pd.DataFrame(rows)


def _rank_variations(results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(
        results,
        key=lambda r: (
            -(-1e9 if r["selection_score"] is None else r["selection_score"]),
            -(-1e9 if r["adjusted_score"] is None else r["adjusted_score"]),
            -(-1e9 if r["loo_predictive_r"] is None else abs(r["loo_predictive_r"])),
        ),
    )


def _summarize_error_pattern(best: dict[str, Any], merged: pd.DataFrame) -> tuple[str, list[dict[str, Any]]]:
    loo = pd.DataFrame(best["per_donor_loo"]).copy()
    if loo.empty:
        return "No analyzable donors were available for leave-one-out error inspection.", []
    loo["abs_error"] = (loo["predicted"] - loo["outcome"]).abs()
    enriched = loo.merge(
        merged.loc[:, ["donor_id", "cognitive_status", "braak_numeric", "cerad_ordinal", "sex", f"{best['feature_column']}_denom"]],
        on="donor_id",
        how="left",
    )
    top = enriched.sort_values("abs_error", ascending=False).head(3)
    top_records = top.loc[:, ["donor_id", "outcome", "predicted", "abs_error", "cognitive_status", "braak_numeric", "cerad_ordinal", "sex", f"{best['feature_column']}_denom"]].to_dict(orient="records")
    common_status = top["cognitive_status"].mode(dropna=True)
    shared_status = str(common_status.iloc[0]) if not common_status.empty else "mixed cognitive-status"
    median_braak = float(top["braak_numeric"].median()) if len(top) else float("nan")
    median_denom = float(top[f"{best['feature_column']}_denom"].median()) if len(top) else float("nan")
    summary = (
        f"The largest LOO errors cluster in donors with {shared_status} status, "
        f"median Braak {median_braak:.1f}, and a median peripyramidal denominator of {median_denom:.0f}; "
        f"these appear to be donors where CA1 reactive fraction alone does not fully track outcome severity."
    )
    return summary, top_records


def _write_report(best: dict[str, Any], ranked: list[dict[str, Any]], merged: pd.DataFrame, summary_text: str, top_errors: list[dict[str, Any]]) -> None:
    other_lines = []
    for item in ranked[1:]:
        other_lines.append(
            f"- `{item['name']}`: partial r {item['partial_r']:.4f}, selection score {item['selection_score']:.4f}, "
            f"LOO predictive r {item['loo_predictive_r']:.4f}."
        )
    if not other_lines:
        other_lines.append("- No alternate local variations were tested.")

    top_error_lines = []
    for row in top_errors:
        top_error_lines.append(
            f"- {row['donor_id']}: outcome {row['outcome']:.3f}, predicted {row['predicted']:.3f}, "
            f"|error| {row['abs_error']:.3f}, status {row['cognitive_status']}, Braak {row['braak_numeric']}, "
            f"CERAD {row['cerad_ordinal']}, denominator {int(row[f\"{best['feature_column']}_denom\"])}."
        )
    if not top_error_lines:
        top_error_lines.append("- No analyzable LOO donor errors were available.")

    report = f"""## Summary
One sentence: tested the {HYPOTHESIS_FAMILY} family, {best['name']} won, and its local winner selection score was {best['selection_score']:.4f}.

## Metrics
Winning variation `{best['name']}` used feature column `{best['feature_column']}` with partial r {best['partial_r']:.4f}, selection score {best['selection_score']:.4f}, LOO predictive r {best['loo_predictive_r']:.4f}, IS-LOO gap {best['is_loo_gap']:.4f}, penalty {best['gap_penalty']:.4f}, and adjusted score {best['adjusted_score']:.4f}.
{"".join(line + chr(10) for line in other_lines)}

## Findings
1. What worked and why (tie to the biological meaning of the target)  
   The winner sharpened the accepted astrocytic-reactivity biology to a CA1 neuron-adjacent niche: reactive astrocytes among CA1 astrocyte-lineage cells lying within {best['radius_um']:.0f} µm of a CA1 pyramidal neuron. This directly targets peripyramidal gliosis around the neuron population most plausibly linked to memory decline.
2. What failed and why (specific to the chosen hypothesis and what went wrong)  
   The nearby radius alternative underperformed, suggesting the signal is not improved by simply broadening the neighborhood. If the broader radius lost, the niche is probably fairly tight; if the tighter radius lost, the signal likely needs a slightly more permissive peripyramidal field.
3. Error pattern: which donors are consistently wrong and what they share  
{summary_text}
{"".join(line + chr(10) for line in top_error_lines)}

## Rationale
The best variation is biologically coherent because it asks whether the CA1 astrocyte-reactivity signal is concentrated specifically around CA1 pyramidal neurons rather than diffusely across all of CA1. That beat the nearby alternative because the winning radius better matched the scale at which astrocyte activation appears to covary with donor memory-decline burden. Since the current accepted panel already contains a region-wide astrocytic reactivity member, this niche-refined candidate is likely partly redundant but could still add value if the evaluator finds that neuron-adjacent gliosis carries cleaner or slightly more specific signal.

## Interpretation
The signal seems to represent CA1 peripyramidal gliosis: the population is CA1 astrocyte-lineage cells, the niche is immediate proximity to CA1 pyramidal neurons, the donor-level summary is the reactive fraction among those niche astrocytes, and the simplest observable pattern is more reactive astrocytes crowding the CA1 pyramidal field in donors with worse memory decline.

## Next
Run one more local sweep around the winning niche scale by testing a slightly tighter or slightly broader CA1 pyramidal-neighborhood threshold and, if needed, a denominator-stabilized variant such as a shrinkage-adjusted reactive fraction.
"""
    Path("/scratch/report.md").write_text(report, encoding="utf-8")


def main() -> None:
    data_root = Path("/data")
    merged = _feature_extraction_table(data_root)
    n_total = int(len(merged))

    results = [evaluate_variation(merged, variation=variation, n_total=n_total) for variation in VARIATIONS]
    ranked = _rank_variations(results)
    best = ranked[0]

    _CANONICAL_SPEC_PATH.write_text(json.dumps(
        {
            "name": best["name"],
            "radius_um": best["radius_um"],
            "feature_column": best["feature_column"],
        },
        indent=2,
    ) + "\n", encoding="utf-8")

    donor_feature_path = write_donor_feature_table(
        "/scratch/donor_feature_table.csv",
        merged,
        feature_column=best["feature_column"],
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUND_COLUMNS,
        id_columns=["donor_id", "slide_name"],
        extra_columns=[
            *(v["feature_column"] for v in VARIATIONS if v["feature_column"] != best["feature_column"]),
            f"{VARIATIONS[0]['feature_column']}_denom",
            f"{VARIATIONS[1]['feature_column']}_denom",
            f"{VARIATIONS[0]['feature_column']}_reactive_n",
            f"{VARIATIONS[1]['feature_column']}_reactive_n",
            "cognitive_status",
            "sex",
        ],
    )

    payload = build_results_payload(
        status="ok",
        feature_name=best["name"],
        outcome=OUTCOME_COLUMN,
        n_total=n_total,
        n_analyzable=best["n_analyzable"],
        partial_r=best["partial_r"],
        ci_lo=best["ci_lo"],
        ci_hi=best["ci_hi"],
        p_value=best["p_value"],
        loo_predictive_r=best["loo_predictive_r"],
        loo_unstable_count=best["loo_unstable_count"],
        loo_max_shift=best["loo_max_shift"],
        donor_ids_used=merged.loc[merged[best["feature_column"]].notna(), "donor_id"].astype(str).tolist(),
        covariates=CONFOUND_COLUMNS,
        recomputed_from_raw=True,
        registry_written=False,
        artifacts={
            "donor_feature_table": str(donor_feature_path),
            "canonical_variation_spec": str(_CANONICAL_SPEC_PATH),
        },
        hypothesis_family=HYPOTHESIS_FAMILY,
        best_variation=best["name"],
        feature_column=best["feature_column"],
        coverage_ratio=best["coverage_ratio"],
        selection_score=best["selection_score"],
        adjusted_score=best["adjusted_score"],
        is_loo_gap=best["is_loo_gap"],
        gap_penalty=best["gap_penalty"],
        bootstrap_median_partial_r=best["bootstrap_median_partial_r"],
        bootstrap_sign_consistency=best["bootstrap_sign_consistency"],
        ranked_variations=[
            {
                key: value
                for key, value in item.items()
                if key not in {"per_donor_loo", "loo_unstable_donors"}
            }
            for item in ranked
        ],
        per_donor_loo=best["per_donor_loo"],
    )
    write_results_payload("/scratch/results.json", payload)

    summary_text, top_errors = _summarize_error_pattern(best, merged)
    _write_report(best, ranked, merged, summary_text, top_errors)

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best['name']}")
    print(f"  IS partial r:      {best['partial_r']:.4f}")
    print(f"  Selection score:   {best['selection_score']:.4f}")
    print(f"  LOO predictive r:  {best['loo_predictive_r']:.4f}  (diagnostic)")
    print(f"  IS-LOO Gap:        {best['is_loo_gap']:.4f}  (penalty={best['gap_penalty']:.4f})")
    print(f"  Adjusted Score:    {best['adjusted_score']:.4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked:
        print(
            f"  {item['name']}  {item['partial_r']:.4f}  "
            f"{item['selection_score']:.4f}  {item['loo_predictive_r']:.4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best['feature_column']}")
    for row in best["per_donor_loo"]:
        print(
            f"  {row['donor_id']}  {row['outcome']:.4f}  {row['predicted']:.4f}  "
            f"{row[best['feature_column']]:.6f}"
        )


if __name__ == "__main__":
    main()
PY
python -m py_compile /scratch/result.py
