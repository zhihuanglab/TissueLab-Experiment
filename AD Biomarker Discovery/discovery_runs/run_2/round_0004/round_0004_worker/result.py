from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, "/shared/lib")

from shared_analysis import build_cell_table, load_training_cohort
from shared_analysis.artifacts import write_donor_feature_table
from shared_analysis.stats import leave_one_out_summary, partial_correlation

DATA_ROOT = Path("/data")
SCRATCH_ROOT = Path("/scratch")
CONTEXT_PATH = SCRATCH_ROOT / "context_bundle.json"
RESULTS_PATH = SCRATCH_ROOT / "results.json"
DONOR_TABLE_PATH = SCRATCH_ROOT / "donor_feature_table.csv"

HYPOTHESIS_FAMILY = "ca1_reactive_astrocyte_hypertrophy"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUNDS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
MIN_CELL_COUNT = 20

VARIATION_SPECS: dict[str, dict[str, str]] = {
    "candidate_variant_a": {
        "description": (
            "Baseline CA1 reactive astrocyte hypertrophy contrast: "
            "median(log1p(area_px2)) of CA1 Reactive Astrocytes minus "
            "median(log1p(area_px2)) of CA1 Astrocytes."
        ),
        "feature_column": "reactive_minus_astrocyte_median_log1p_area_ca1",
        "reactive_summary_column": "reactive_median_log1p_area_ca1",
        "astrocyte_summary_column": "astrocyte_median_log1p_area_ca1",
    },
    "candidate_variant_b": {
        "description": (
            "Upper-tail CA1 reactive astrocyte hypertrophy contrast: "
            "75th percentile(log1p(area_px2)) of CA1 Reactive Astrocytes minus "
            "75th percentile(log1p(area_px2)) of CA1 Astrocytes."
        ),
        "feature_column": "reactive_minus_astrocyte_p75_log1p_area_ca1",
        "reactive_summary_column": "reactive_p75_log1p_area_ca1",
        "astrocyte_summary_column": "astrocyte_p75_log1p_area_ca1",
    },
}

# This value is rewritten after the sweep finishes so held-out replay uses the winner.
CANONICAL_REPLAY_VARIATION = "candidate_variant_a"


def _canonical_variation_name() -> str:
    if CANONICAL_REPLAY_VARIATION in VARIATION_SPECS:
        return CANONICAL_REPLAY_VARIATION
    return "candidate_variant_a"


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


def _design_matrix(frame: pd.DataFrame) -> np.ndarray:
    arr = frame.to_numpy(dtype=float)
    intercept = np.ones((len(frame), 1), dtype=float)
    return np.concatenate([intercept, arr], axis=1)


def residualized_loo_predictions(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
    id_col: str = "donor_id",
) -> tuple[float, pd.DataFrame]:
    cols = [id_col, feature_col, outcome_col, *confounds]
    frame = df.loc[:, cols].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    predictions: list[dict[str, float | str]] = []
    pred_values: list[float] = []
    actual_values: list[float] = []
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
        if denom <= 0:
            pred = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            resid_feature_test = test[feature_col].to_numpy(dtype=float) - x_test_conf @ beta_feature
            pred = float(slope * resid_feature_test[0])

        actual = float(
            test[outcome_col].to_numpy(dtype=float)[0] - (x_test_conf @ beta_outcome)[0]
        )
        pred_values.append(pred)
        actual_values.append(actual)
        predictions.append(
            {
                id_col: str(test.iloc[0][id_col]),
                "outcome": float(test.iloc[0][outcome_col]),
                "predicted": pred,
                "actual_residualized_outcome": actual,
            }
        )

    pred_arr = np.asarray(pred_values, dtype=float)
    act_arr = np.asarray(actual_values, dtype=float)
    mask = np.isfinite(pred_arr) & np.isfinite(act_arr)
    loo_r = float("nan")
    if mask.sum() >= 3 and np.std(pred_arr[mask]) > 0 and np.std(act_arr[mask]) > 0:
        loo_r = float(np.corrcoef(pred_arr[mask], act_arr[mask])[0, 1])

    return loo_r, pd.DataFrame(predictions)


def _gap_penalty(is_loo_gap: float | None) -> tuple[float | None, float | None]:
    gap = _safe_float(is_loo_gap)
    if gap is None:
        return None, None
    penalty = max(0.0, gap - 0.15) * 0.5
    if gap > 0.30:
        adjusted = -1.0
    else:
        adjusted = None
    return penalty, adjusted


def extract_donor_features_from_slide(slide_path: str | Path) -> dict[str, float | int | None]:
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    slide_path = Path(slide_path)
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=True)
    ca1 = cells.loc[
        (cells["region"] == "CA1")
        & (cells["cell_type"].isin(["Reactive Astrocyte", "Astrocyte"])),
        ["cell_type", "area"],
    ].copy()

    ca1["area"] = pd.to_numeric(ca1["area"], errors="coerce")
    ca1 = ca1.loc[np.isfinite(ca1["area"])].copy()
    ca1["log1p_area"] = np.log1p(np.clip(ca1["area"].to_numpy(dtype=float), a_min=0.0, a_max=None))

    reactive = ca1.loc[ca1["cell_type"] == "Reactive Astrocyte", "log1p_area"].to_numpy(dtype=float)
    astrocyte = ca1.loc[ca1["cell_type"] == "Astrocyte", "log1p_area"].to_numpy(dtype=float)

    result: dict[str, float | int | None] = {
        "ca1_reactive_astrocyte_count": int(len(reactive)),
        "ca1_astrocyte_count": int(len(astrocyte)),
        "ca1_reactive_astrocyte_plus_astrocyte_count": int(len(ca1)),
    }

    if len(reactive) < MIN_CELL_COUNT or len(astrocyte) < MIN_CELL_COUNT:
        for spec in VARIATION_SPECS.values():
            result[spec["reactive_summary_column"]] = None
            result[spec["astrocyte_summary_column"]] = None
            result[spec["feature_column"]] = None
        return result

    reactive_median = float(np.median(reactive))
    astro_median = float(np.median(astrocyte))
    reactive_p75 = float(np.quantile(reactive, 0.75))
    astro_p75 = float(np.quantile(astrocyte, 0.75))

    result.update(
        {
            "reactive_median_log1p_area_ca1": reactive_median,
            "astrocyte_median_log1p_area_ca1": astro_median,
            "reactive_minus_astrocyte_median_log1p_area_ca1": reactive_median - astro_median,
            "reactive_p75_log1p_area_ca1": reactive_p75,
            "astrocyte_p75_log1p_area_ca1": astro_p75,
            "reactive_minus_astrocyte_p75_log1p_area_ca1": reactive_p75 - astro_p75,
        }
    )
    return result


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    features = extract_donor_features_from_slide(data_root / slide_name)
    feature_col = VARIATION_SPECS[_canonical_variation_name()]["feature_column"]
    value = features.get(feature_col)
    return None if value is None else float(value)


def build_feature_table(data_root: str | Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = (cohort["sex"].astype(str).str.lower() == "male").astype(float)

    feature_rows = []
    for _, row in cohort.iterrows():
        feats = extract_donor_features_from_slide(Path(data_root) / str(row["slide_name"]))
        feats["donor_id"] = str(row["donor_id"])
        feature_rows.append(feats)
    feature_df = pd.DataFrame(feature_rows)

    merged = cohort.merge(feature_df, on="donor_id", how="left")
    return merged


def evaluate_variation(
    table: pd.DataFrame,
    *,
    variation_name: str,
) -> dict[str, Any]:
    spec = VARIATION_SPECS[variation_name]
    feature_col = spec["feature_column"]
    partial = partial_correlation(
        table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
    )
    loo_r, loo_df = residualized_loo_predictions(
        table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
        id_col="donor_id",
    )
    full_df = table.loc[
        :,
        ["donor_id", "slide_name", OUTCOME_COLUMN, feature_col, *CONFOUNDS],
    ].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    merge_columns = [
        "donor_id",
        feature_col,
        spec["reactive_summary_column"],
        spec["astrocyte_summary_column"],
        "ca1_reactive_astrocyte_count",
        "ca1_astrocyte_count",
    ]
    for optional_col in ["sex", "cognitive_status", "overall_ad_neuropath_change"]:
        if optional_col in table.columns:
            merge_columns.append(optional_col)
    loo_df = loo_df.merge(
        table.loc[:, merge_columns],
        on="donor_id",
        how="left",
    )

    n_total = int(len(table))
    n_analyzable = int(partial["n"]) if _safe_float(partial.get("n")) is not None else 0
    partial_r = _safe_float(partial.get("partial_r"))
    selection_score = None
    if partial_r is not None and n_total > 0:
        selection_score = abs(partial_r) * (n_analyzable / n_total)

    is_loo_gap = None
    if partial_r is not None and _safe_float(loo_r) is not None:
        is_loo_gap = abs(partial_r) - float(loo_r)
    penalty, forced_adjusted = _gap_penalty(is_loo_gap)
    if forced_adjusted is not None:
        adjusted_score = forced_adjusted
    else:
        adjusted_score = None if _safe_float(loo_r) is None else float(loo_r) - (penalty or 0.0)

    instability = leave_one_out_summary(
        full_df,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
        id_col="donor_id",
    )

    return {
        "variation_name": variation_name,
        "description": spec["description"],
        "feature_column": feature_col,
        "feature_name": f"{HYPOTHESIS_FAMILY}__{feature_col}",
        "reactive_summary_column": spec["reactive_summary_column"],
        "astrocyte_summary_column": spec["astrocyte_summary_column"],
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "partial_r": partial_r,
        "p_value": _safe_float(partial.get("p_value")),
        "selection_score": selection_score,
        "loo_predictive_r": _safe_float(loo_r),
        "is_loo_gap": is_loo_gap,
        "penalty": penalty,
        "adjusted_score": adjusted_score,
        "loo_unstable_count": int(instability.get("unstable_count", 0) or 0),
        "loo_max_shift": _safe_float(instability.get("max_shift")),
        "per_donor_loo": loo_df.to_dict(orient="records"),
    }


def _fmt(value: Any) -> str:
    val = _safe_float(value)
    return "nan" if val is None else f"{val:.4f}"


def main() -> None:
    warnings.filterwarnings("ignore")
    context = json.loads(CONTEXT_PATH.read_text())
    brief = context["worker_brief"]
    planned_variations = [row["name"] for row in brief["variations"]]
    table = build_feature_table(DATA_ROOT)

    metrics = [evaluate_variation(table, variation_name=name) for name in planned_variations]
    ranked = sorted(
        metrics,
        key=lambda row: (
            -999.0 if row["selection_score"] is None else row["selection_score"],
            -999.0 if row["adjusted_score"] is None else row["adjusted_score"],
        ),
        reverse=True,
    )
    best = ranked[0]

    extra_cols = list(
        dict.fromkeys(
            [
                *(spec["feature_column"] for spec in VARIATION_SPECS.values()),
                *(spec["reactive_summary_column"] for spec in VARIATION_SPECS.values()),
                *(spec["astrocyte_summary_column"] for spec in VARIATION_SPECS.values()),
                "ca1_reactive_astrocyte_count",
                "ca1_astrocyte_count",
                "ca1_reactive_astrocyte_plus_astrocyte_count",
                "sex",
            ]
        )
    )
    write_donor_feature_table(
        DONOR_TABLE_PATH,
        table,
        feature_column=best["feature_column"],
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUNDS,
        extra_columns=extra_cols,
    )

    results_payload = {
        "status": "ok",
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "scientific_question": brief["scientific_question"],
        "approach": brief["approach"],
        "best_variation": best["variation_name"],
        "feature_name": best["feature_name"],
        "feature_column": best["feature_column"],
        "n_total": best["n_total"],
        "n_analyzable": best["n_analyzable"],
        "partial_r": best["partial_r"],
        "selection_score": best["selection_score"],
        "loo_predictive_r": best["loo_predictive_r"],
        "is_loo_gap": best["is_loo_gap"],
        "penalty": best["penalty"],
        "adjusted_score": best["adjusted_score"],
        "p_value": best["p_value"],
        "loo_unstable_count": best["loo_unstable_count"],
        "loo_max_shift": best["loo_max_shift"],
        "ranked_variations": ranked,
        "canonical_replay_variation": best["variation_name"],
        "artifacts": {
            "donor_feature_table": str(DONOR_TABLE_PATH),
        },
    }
    RESULTS_PATH.write_text(json.dumps(results_payload, indent=2) + "\n", encoding="utf-8")

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best['variation_name']}")
    print(f"  IS partial r:      {_fmt(best['partial_r'])}")
    print(f"  Selection score:   {_fmt(best['selection_score'])}")
    print(f"  LOO predictive r:  {_fmt(best['loo_predictive_r'])}  (diagnostic)")
    print(
        f"  IS-LOO Gap:        {_fmt(best['is_loo_gap'])}  "
        f"(penalty={_fmt(best['penalty'])})"
    )
    print(f"  Adjusted Score:    {_fmt(best['adjusted_score'])}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for row in ranked:
        print(
            f"  {row['variation_name']}  {_fmt(row['partial_r'])}  "
            f"{_fmt(row['selection_score'])}  {_fmt(row['loo_predictive_r'])}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(
        "  donor_id  outcome  predicted  "
        f"{best['feature_column']}  {best['reactive_summary_column']}  "
        f"{best['astrocyte_summary_column']}  ca1_reactive_astrocyte_count  ca1_astrocyte_count"
    )
    for row in best["per_donor_loo"]:
        print(
            f"  {row['donor_id']}  {_fmt(row['outcome'])}  {_fmt(row['predicted'])}  "
            f"{_fmt(row.get(best['feature_column']))}  "
            f"{_fmt(row.get(best['reactive_summary_column']))}  "
            f"{_fmt(row.get(best['astrocyte_summary_column']))}  "
            f"{int(row.get('ca1_reactive_astrocyte_count', 0))}  "
            f"{int(row.get('ca1_astrocyte_count', 0))}"
        )


if __name__ == "__main__":
    main()
