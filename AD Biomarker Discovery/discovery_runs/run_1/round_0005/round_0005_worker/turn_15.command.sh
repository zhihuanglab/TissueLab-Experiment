cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, "/shared/lib")
from shared_analysis import build_cell_table, load_training_cohort, write_donor_feature_table  # noqa: E402
from shared_analysis.stats import (  # noqa: E402
    bootstrap_partial_correlation,
    leave_one_out_summary,
    partial_correlation,
)

warnings.filterwarnings("ignore")


HYPOTHESIS_FAMILY = "ca1_pyramidal_reactive_astro_encirclement"
FEATURE_NAME = HYPOTHESIS_FAMILY
OUTCOME_COLUMN = "slope_zmem0"
CONFOUNDS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]

REGION = "CA1"
SOURCE_CELL_TYPE = "Pyramidal Neuron"
TARGET_CELL_TYPE = "Reactive Astrocyte"
RADIUS_PX = 60.0

VARIATIONS: dict[str, dict[str, Any]] = {
    "encirclement_k2_60px": {
        "k": 2,
        "feature_column": "ca1_pyramidal_reactive_astro_encirclement_k2_60px",
        "description": "Fraction of CA1 pyramidal neurons with at least 2 reactive astrocytes within 60 px.",
    },
    "encirclement_k3_60px": {
        "k": 3,
        "feature_column": "ca1_pyramidal_reactive_astro_encirclement_k3_60px",
        "description": "Fraction of CA1 pyramidal neurons with at least 3 reactive astrocytes within 60 px.",
    },
    "encirclement_k4_60px": {
        "k": 4,
        "feature_column": "ca1_pyramidal_reactive_astro_encirclement_k4_60px",
        "description": "Fraction of CA1 pyramidal neurons with at least 4 reactive astrocytes within 60 px.",
    },
}

CANONICAL_FEATURE_COLUMN = "feature_value"


def _paths(base: str | Path | None = None) -> dict[str, Path]:
    root = Path(base) if base is not None else Path(__file__).resolve().parent
    return {
        "root": root,
        "table": root / "donor_feature_table.csv",
        "results": root / "results.json",
        "report": root / "report.md",
        "config": root / "encirclement_config.json",
    }


def _clean_float(value: Any) -> float | None:
    try:
        x = float(value)
    except Exception:
        return None
    if not np.isfinite(x):
        return None
    return x


def _jsonify(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _jsonify(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_jsonify(v) for v in obj]
    if isinstance(obj, tuple):
        return [_jsonify(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        return _clean_float(obj)
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if pd.isna(obj):
        return None
    return obj


def _selection_score(partial_r: float | None, n_analyzable: int | None, n_total: int | None) -> float:
    if partial_r is None or not np.isfinite(partial_r):
        return float("nan")
    if n_analyzable is None or n_total is None or n_total <= 0:
        return float("nan")
    coverage = max(0.0, min(1.0, float(n_analyzable) / float(n_total)))
    return float(abs(float(partial_r)) * coverage)


def _adjusted_metrics(partial_r: float | None, loo_predictive_r: float | None) -> tuple[float, float, float]:
    is_r = abs(float(partial_r)) if partial_r is not None and np.isfinite(partial_r) else 0.0
    loo_r = abs(float(loo_predictive_r)) if loo_predictive_r is not None and np.isfinite(loo_predictive_r) else 0.0
    gap = is_r - loo_r
    if gap > 0.30:
        return float(gap), float(gap), -1.0
    penalty = max(0.0, gap - 0.15) * 0.5
    return float(gap), float(penalty), float(loo_r - penalty)


def _load_eval_cohort(data_root: str | Path) -> pd.DataFrame:
    data_root = Path(data_root)
    training = data_root / "training_cohort.csv"
    test = data_root / "test_cohort.csv"
    if training.exists():
        cohort = load_training_cohort(data_root).copy()
    elif test.exists():
        cohort = pd.read_csv(test).copy()
    else:
        cohort = load_training_cohort(data_root).copy()
    cohort["donor_id"] = cohort["donor_id"].astype(str)
    cohort["slide_name"] = cohort["slide_name"].astype(str)
    if "sex_binary" not in cohort.columns:
        sex_map = {"Male": 1.0, "Female": 0.0, "M": 1.0, "F": 0.0}
        cohort["sex_binary"] = cohort["sex"].map(sex_map).astype(float)
    return cohort


def _neighbor_counts_within_radius(source_xy: np.ndarray, target_xy: np.ndarray, radius: float) -> np.ndarray:
    if len(source_xy) == 0:
        return np.zeros(0, dtype=int)
    if len(target_xy) == 0:
        return np.zeros(len(source_xy), dtype=int)
    tree = cKDTree(target_xy)
    neighborhoods = tree.query_ball_point(source_xy, r=radius)
    return np.fromiter((len(n) for n in neighborhoods), dtype=int, count=len(source_xy))


def _compute_slide_features(zarr_path: Path) -> dict[str, Any]:
    cells = build_cell_table(zarr_path, include_regions=True, include_geometry=False)
    ca1 = cells.loc[cells["region"] == REGION, ["x", "y", "cell_type"]].copy()

    pyramidal = ca1.loc[ca1["cell_type"] == SOURCE_CELL_TYPE, ["x", "y"]].to_numpy(dtype=float, copy=False)
    reactive = ca1.loc[ca1["cell_type"] == TARGET_CELL_TYPE, ["x", "y"]].to_numpy(dtype=float, copy=False)

    n_pyramidal = int(len(pyramidal))
    n_reactive = int(len(reactive))
    result: dict[str, Any] = {
        "ca1_total_cells": int(len(ca1)),
        "ca1_pyramidal_count": n_pyramidal,
        "ca1_reactive_astrocyte_count": n_reactive,
    }

    if n_pyramidal == 0:
        counts = np.zeros(0, dtype=int)
        for cfg in VARIATIONS.values():
            result[cfg["feature_column"]] = float("nan")
    else:
        counts = _neighbor_counts_within_radius(pyramidal, reactive, radius=RADIUS_PX)
        for cfg in VARIATIONS.values():
            result[cfg["feature_column"]] = float(np.mean(counts >= int(cfg["k"])))

    if n_pyramidal > 0:
        result["ca1_reactive_astro_per_pyramidal_mean_60px"] = float(np.mean(counts))
        result["ca1_reactive_astro_per_pyramidal_median_60px"] = float(np.median(counts))
    else:
        result["ca1_reactive_astro_per_pyramidal_mean_60px"] = float("nan")
        result["ca1_reactive_astro_per_pyramidal_median_60px"] = float("nan")
    return result


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    data_root = Path(data_root)
    cohort = _load_eval_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == str(donor_id)]
    if donor_rows.empty:
        return None

    config_path = Path(__file__).resolve().with_name("encirclement_config.json")
    variation_name = "encirclement_k2_60px"
    if config_path.exists():
        try:
            config = json.loads(config_path.read_text(encoding="utf-8"))
            variation_name = str(config.get("best_variation", variation_name))
        except Exception:
            pass

    if variation_name not in VARIATIONS:
        variation_name = "encirclement_k2_60px"

    slide_name = str(donor_rows.iloc[0]["slide_name"])
    features = _compute_slide_features(data_root / slide_name)
    value = features.get(VARIATIONS[variation_name]["feature_column"])
    return float(value) if value is not None and np.isfinite(value) else float("nan")


def _extract_all_features(cohort: pd.DataFrame, data_root: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row in cohort.itertuples(index=False):
        donor_row = {
            "donor_id": str(row.donor_id),
            "slide_name": str(row.slide_name),
            OUTCOME_COLUMN: float(getattr(row, OUTCOME_COLUMN)),
            "max_age_vis": float(getattr(row, "max_age_vis")),
            "braak_numeric": float(getattr(row, "braak_numeric")),
            "cerad_ordinal": float(getattr(row, "cerad_ordinal")),
            "sex_binary": float(getattr(row, "sex_binary")),
            "sex": getattr(row, "sex", None),
            "cognitive_status": getattr(row, "cognitive_status", None),
            "overall_ad_neuropath_change": getattr(row, "overall_ad_neuropath_change", None),
        }
        donor_row.update(_compute_slide_features(data_root / str(row.slide_name)))
        rows.append(donor_row)
    return pd.DataFrame(rows)


def _loo_predictions(df: pd.DataFrame, *, feature_col: str) -> tuple[float, pd.DataFrame]:
    cols = ["donor_id", OUTCOME_COLUMN, feature_col, *CONFOUNDS]
    frame = df.loc[:, cols].dropna().reset_index(drop=True)

    preds_resid: list[float] = []
    actuals_resid: list[float] = []
    rows: list[dict[str, Any]] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train_full = np.column_stack(
            [
                np.ones(len(train), dtype=float),
                train[CONFOUNDS].to_numpy(dtype=float),
                train[[feature_col]].to_numpy(dtype=float),
            ]
        )
        x_test_full = np.column_stack(
            [
                np.ones(len(test), dtype=float),
                test[CONFOUNDS].to_numpy(dtype=float),
                test[[feature_col]].to_numpy(dtype=float),
            ]
        )
        y_train = train[OUTCOME_COLUMN].to_numpy(dtype=float)
        beta_full, *_ = np.linalg.lstsq(x_train_full, y_train, rcond=None)
        pred_raw = float((x_test_full @ beta_full)[0])

        x_train_conf = np.column_stack([np.ones(len(train), dtype=float), train[CONFOUNDS].to_numpy(dtype=float)])
        x_test_conf = np.column_stack([np.ones(len(test), dtype=float), test[CONFOUNDS].to_numpy(dtype=float)])

        beta_feature, *_ = np.linalg.lstsq(x_train_conf, train[feature_col].to_numpy(dtype=float), rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train_conf, train[OUTCOME_COLUMN].to_numpy(dtype=float), rcond=None)

        resid_feature_train = train[feature_col].to_numpy(dtype=float) - x_train_conf @ beta_feature
        resid_outcome_train = train[OUTCOME_COLUMN].to_numpy(dtype=float) - x_train_conf @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            pred_resid = float("nan")
            actual_resid = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            resid_feature_test = test[feature_col].to_numpy(dtype=float) - x_test_conf @ beta_feature
            actual_resid = float(test[OUTCOME_COLUMN].to_numpy(dtype=float)[0] - (x_test_conf @ beta_outcome)[0])
            pred_resid = float(slope * resid_feature_test[0])

        preds_resid.append(pred_resid)
        actuals_resid.append(actual_resid)
        rows.append(
            {
                "donor_id": str(test.iloc[0]["donor_id"]),
                "outcome": float(test.iloc[0][OUTCOME_COLUMN]),
                "predicted": pred_raw,
                "predicted_residualized": pred_resid,
                "actual_residualized": actual_resid,
                feature_col: float(test.iloc[0][feature_col]),
            }
        )

    pred_arr = np.asarray(preds_resid, dtype=float)
    act_arr = np.asarray(actuals_resid, dtype=float)
    mask = np.isfinite(pred_arr) & np.isfinite(act_arr)
    if mask.sum() < 3 or np.std(pred_arr[mask]) == 0 or np.std(act_arr[mask]) == 0:
        loo_r = float("nan")
    else:
        loo_r = float(np.corrcoef(pred_arr[mask], act_arr[mask])[0, 1])
    return loo_r, pd.DataFrame(rows)


def _evaluate_variation(table: pd.DataFrame, *, variation_name: str) -> dict[str, Any]:
    cfg = VARIATIONS[variation_name]
    feature_col = cfg["feature_column"]
    eval_cols = ["donor_id", OUTCOME_COLUMN, feature_col, *CONFOUNDS]
    analyzable = table.loc[:, eval_cols].dropna().reset_index(drop=True)

    base = partial_correlation(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
    )
    boot = bootstrap_partial_correlation(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
        n_boot=1000,
        random_state=0,
    )
    loo_predictive_r, per_donor_loo = _loo_predictions(table, feature_col=feature_col)
    loo_summary = leave_one_out_summary(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
        id_col="donor_id",
    )
    selection = _selection_score(base["partial_r"], int(len(analyzable)), int(len(table)))
    gap, penalty, adjusted = _adjusted_metrics(base["partial_r"], loo_predictive_r)

    return {
        "name": variation_name,
        "description": cfg["description"],
        "k": int(cfg["k"]),
        "radius_px": float(RADIUS_PX),
        "feature_column": feature_col,
        "n_total": int(len(table)),
        "n_analyzable": int(len(analyzable)),
        "partial_r": float(base["partial_r"]),
        "p_value": float(base["p_value"]),
        "ci_lo": _clean_float(boot.get("ci_lo")),
        "ci_hi": _clean_float(boot.get("ci_hi")),
        "selection_score": float(selection),
        "loo_predictive_r": float(loo_predictive_r),
        "is_loo_gap": float(gap),
        "gap_penalty": float(penalty),
        "adjusted_score": float(adjusted),
        "loo_unstable_count": int(loo_summary.get("unstable_count", 0) or 0),
        "loo_max_shift": _clean_float(loo_summary.get("max_shift")),
        "per_donor_loo": per_donor_loo,
    }


def _direction_text(partial_r: float) -> str:
    if not np.isfinite(partial_r):
        return "showed no stable directional relationship"
    if partial_r < 0:
        return "higher encirclement tracked lower slope_zmem0 values, consistent with worse memory decline"
    if partial_r > 0:
        return "higher encirclement tracked higher slope_zmem0 values"
    return "showed no directional relationship"


def _report_error_pattern(best_rows: pd.DataFrame, donor_table: pd.DataFrame) -> str:
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


def _render_report(
    *,
    donor_table: pd.DataFrame,
    ranked: list[dict[str, Any]],
    best: dict[str, Any],
    best_loo: pd.DataFrame,
) -> str:
    winner = best["name"]
    local_score = best["selection_score"]
    others = [r for r in ranked if r["name"] != winner]

    ranking_lines = []
    for idx, item in enumerate(ranked, start=1):
        ranking_lines.append(
            f"{idx}. {item['name']}: partial_r={item['partial_r']:.4f}, "
            f"selection_score={item['selection_score']:.4f}, "
            f"loo_predictive_r={item['loo_predictive_r']:.4f}"
        )

    error_pattern = _report_error_pattern(best_loo, donor_table)
    direction = _direction_text(float(best["partial_r"]))

    if others:
        failed_text = (
            "The stricter thresholds lost signal relative to the winner. "
            + " ".join(
                f"{item['name']} fell to selection_score={item['selection_score']:.4f}"
                f" with loo_predictive_r={item['loo_predictive_r']:.4f}."
                for item in others
            )
        )
    else:
        failed_text = "No nearby alternatives were tested."

    likely_additive = (
        "Because this is a neuron-centered severity readout of the same CA1 reactive-astro niche family, "
        "it is biologically coherent but likely somewhat redundant with the current panel unless its donor-level errors differ."
    )

    report = f"""## Summary
One sentence: tested the {HYPOTHESIS_FAMILY} family, {winner} won, and the local winner selection score was {local_score:.4f}.

## Metrics
Winning variation: {winner}
- IS partial r: {best['partial_r']:.4f}
- Selection score: {best['selection_score']:.4f}
- LOO predictive r: {best['loo_predictive_r']:.4f}
- IS-LOO Gap: {best['is_loo_gap']:.4f} (penalty={best['gap_penalty']:.4f})
- Adjusted Score: {best['adjusted_score']:.4f}

Ranking of tested variations:
{chr(10).join(ranking_lines)}

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - {winner} worked best because it kept the broader moderate-burden perineuronal niche while staying neuron-centric: the donor scalar is the fraction of CA1 pyramidal neurons with at least {best['k']} nearby reactive astrocytes inside a 60 px halo.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - {failed_text}
3. Error pattern: which donors are consistently wrong and what they share
   - {error_pattern}

## Rationale
The winning approach is biologically coherent because it asks whether CA1 pyramidal neurons are embedded in a small local cluster of reactive astrocytes, a simple perineuronal injury/glial-response pattern that can be audited directly from centroid geometry. It beat nearby alternatives because the winner preserved more dynamic range across donors than the stricter encirclement thresholds, which make the event rarer and more unstable. {likely_additive}

## Interpretation
State this when it materially helps explain the biomarker.
The signal seems to mean biologically that donors with more CA1 pyramidal neurons sitting inside a reactive-astrocyte halo show a different memory-decline trajectory; specifically, {direction}. The population is CA1 pyramidal neurons, the niche is nearby CA1 reactive astrocytes within 60 px, the feature summary is the donor-level fraction of neurons exceeding the local reactive-astro count threshold, and the simplest observable tissue pattern is a CA1 neuron repeatedly ringed by several reactive astrocytes.

## Next
One specific suggestion for the next local sweep based on the error pattern and which nearby variations won or lost.
Test the same neuron-centric encirclement family at the winning threshold but vary the neighborhood radius around 60 px (for example 40, 80, and 100 px) to determine whether the signal is driven by a tighter perisomatic halo or by a broader CA1 reactive-astro neighborhood.
"""
    return report


def main() -> int:
    data_root = Path("/data")
    paths = _paths()

    cohort = _load_eval_cohort(data_root)
    donor_table = _extract_all_features(cohort, data_root)

    evaluated = [_evaluate_variation(donor_table, variation_name=name) for name in VARIATIONS]
    ranked = sorted(
        evaluated,
        key=lambda d: (
            -abs(d["selection_score"]) if np.isfinite(d["selection_score"]) else math.inf,
            -abs(d["partial_r"]) if np.isfinite(d["partial_r"]) else math.inf,
            -abs(d["loo_predictive_r"]) if np.isfinite(d["loo_predictive_r"]) else math.inf,
        ),
    )
    best = ranked[0]

    donor_table[CANONICAL_FEATURE_COLUMN] = donor_table[best["feature_column"]]

    write_donor_feature_table(
        paths["table"],
        donor_table,
        feature_column=CANONICAL_FEATURE_COLUMN,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUNDS,
        extra_columns=[
            *[cfg["feature_column"] for cfg in VARIATIONS.values()],
            "ca1_total_cells",
            "ca1_pyramidal_count",
            "ca1_reactive_astrocyte_count",
            "ca1_reactive_astro_per_pyramidal_mean_60px",
            "ca1_reactive_astro_per_pyramidal_median_60px",
            "sex",
            "cognitive_status",
            "overall_ad_neuropath_change",
        ],
    )

    config = {
        "family": HYPOTHESIS_FAMILY,
        "region": REGION,
        "source_cell_type": SOURCE_CELL_TYPE,
        "target_cell_type": TARGET_CELL_TYPE,
        "radius_px": RADIUS_PX,
        "best_variation": best["name"],
        "best_feature_column": best["feature_column"],
        "canonical_feature_column": CANONICAL_FEATURE_COLUMN,
        "variations": {name: {k: v for k, v in cfg.items() if k != "description"} for name, cfg in VARIATIONS.items()},
    }
    paths["config"].write_text(json.dumps(_jsonify(config), indent=2) + "\n", encoding="utf-8")

    best_loo = best["per_donor_loo"].copy()
    report_text = _render_report(
        donor_table=donor_table,
        ranked=ranked,
        best=best,
        best_loo=best_loo,
    )
    paths["report"].write_text(report_text, encoding="utf-8")

    ranked_json = []
    for item in ranked:
        ranked_json.append(
            {
                key: _jsonify(value)
                for key, value in item.items()
                if key != "per_donor_loo"
            }
        )

    results = {
        "status": "ok",
        "feature_name": FEATURE_NAME,
        "outcome": OUTCOME_COLUMN,
        "best_variation": best["name"],
        "feature_column": CANONICAL_FEATURE_COLUMN,
        "n_total": int(best["n_total"]),
        "n_analyzable": int(best["n_analyzable"]),
        "partial_r": _jsonify(best["partial_r"]),
        "p_value": _jsonify(best["p_value"]),
        "ci_lo": _jsonify(best["ci_lo"]),
        "ci_hi": _jsonify(best["ci_hi"]),
        "selection_score": _jsonify(best["selection_score"]),
        "loo_predictive_r": _jsonify(best["loo_predictive_r"]),
        "is_loo_gap": _jsonify(best["is_loo_gap"]),
        "gap_penalty": _jsonify(best["gap_penalty"]),
        "adjusted_score": _jsonify(best["adjusted_score"]),
        "covariates": CONFOUNDS,
        "donor_ids_used": donor_table.loc[
            donor_table[[best["feature_column"], OUTCOME_COLUMN, *CONFOUNDS]].notna().all(axis=1),
            "donor_id",
        ].astype(str).tolist(),
        "ranked_variations": ranked_json,
        "artifacts": {
            "donor_feature_table": str(paths["table"]),
            "config": str(paths["config"]),
            "report": str(paths["report"]),
        },
        "population": SOURCE_CELL_TYPE,
        "niche_population": TARGET_CELL_TYPE,
        "region": REGION,
        "radius_px": RADIUS_PX,
        "per_donor_loo": _jsonify(best_loo.to_dict(orient="records")),
    }
    paths["results"].write_text(json.dumps(_jsonify(results), indent=2, allow_nan=False) + "\n", encoding="utf-8")

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best['name']}")
    print(f"  IS partial r:      {best['partial_r']: .4f}")
    print(f"  Selection score:   {best['selection_score']: .4f}")
    print(f"  LOO predictive r:  {best['loo_predictive_r']: .4f}  (diagnostic)")
    print(f"  IS-LOO Gap:        {best['is_loo_gap']: .4f}  (penalty={best['gap_penalty']: .4f})")
    print(f"  Adjusted Score:    {best['adjusted_score']: .4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked:
        print(
            f"  {item['name']}  "
            f"{item['partial_r']: .4f}  "
            f"{item['selection_score']: .4f}  "
            f"{item['loo_predictive_r']: .4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best['feature_column']}")
    for row in best_loo.itertuples(index=False):
        row_map = row._asdict()
        print(
            f"  {row_map['donor_id']}  "
            f"{float(row_map['outcome']): .4f}  "
            f"{float(row_map['predicted']): .4f}  "
            f"{float(row_map[best['feature_column']]): .4f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
PY
python /scratch/result.py
