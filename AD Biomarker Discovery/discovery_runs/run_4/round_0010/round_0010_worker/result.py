from __future__ import annotations

import json
import math
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, "/shared/lib")

from shared_analysis import build_cell_table, load_training_cohort
from shared_analysis.artifacts import (
    build_results_payload,
    write_donor_feature_table,
    write_results_payload,
)
from shared_analysis.stats import leave_one_out_summary, partial_correlation


warnings.filterwarnings(
    "ignore",
    message="Object at .* is not recognized as a component of a Zarr hierarchy.",
)

HYPOTHESIS_FAMILY = "ca1-reactive-lymphocyte-niche"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUNDS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
BEST_VARIATION_SIDECAR = Path(__file__).with_name("best_variation.json")

VARIATIONS = [
    {
        "variation_name": "candidate_variant_a",
        "description": (
            "ca1_reactive_lymphocyte_neighbor_fraction_r96px: fraction of CA1 "
            "pyramidal neurons with at least one CA1 reactive astrocyte and at least "
            "one CA1 lymphocyte within 96 px."
        ),
        "feature_column": "ca1_reactive_lymphocyte_neighbor_fraction_r96px",
        "radius_px": 96.0,
    },
    {
        "variation_name": "candidate_variant_b",
        "description": (
            "ca1_reactive_lymphocyte_neighbor_fraction_r128px: same metric with a 128 px "
            "neighborhood to test a slightly broader immune niche."
        ),
        "feature_column": "ca1_reactive_lymphocyte_neighbor_fraction_r128px",
        "radius_px": 128.0,
    },
]
VARIATION_BY_NAME = {v["variation_name"]: v for v in VARIATIONS}
DEFAULT_CANONICAL_VARIATION = VARIATIONS[0]["variation_name"]


@dataclass
class DonorFeatureBundle:
    features: dict[str, float]
    counts: dict[str, int]


def _design_matrix(df: pd.DataFrame) -> np.ndarray:
    x = df.to_numpy(dtype=float)
    intercept = np.ones((len(df), 1), dtype=float)
    return np.concatenate([intercept, x], axis=1)


def _clean_analysis_frame(
    df: pd.DataFrame,
    *,
    feature_col: str,
) -> pd.DataFrame:
    required = ["donor_id", OUTCOME_COLUMN, feature_col, *CONFOUNDS]
    return df.loc[:, required].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)


def _score_to_float(value) -> float:
    try:
        out = float(value)
    except Exception:
        return float("nan")
    return out


def _safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3:
        return float("nan")
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        mask = np.isfinite(x) & np.isfinite(y)
        x = x[mask]
        y = y[mask]
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def _residualized_loo_predictions(
    df: pd.DataFrame,
    *,
    feature_col: str,
) -> tuple[float, pd.DataFrame]:
    frame = _clean_analysis_frame(df, feature_col=feature_col)
    rows: list[dict[str, float | str]] = []
    resid_actuals: list[float] = []
    resid_preds: list[float] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train_conf = _design_matrix(train[CONFOUNDS])
        x_test_conf = _design_matrix(test[CONFOUNDS])

        beta_feature, *_ = np.linalg.lstsq(
            x_train_conf,
            train[feature_col].to_numpy(dtype=float),
            rcond=None,
        )
        beta_outcome, *_ = np.linalg.lstsq(
            x_train_conf,
            train[OUTCOME_COLUMN].to_numpy(dtype=float),
            rcond=None,
        )

        resid_feature_train = train[feature_col].to_numpy(dtype=float) - x_train_conf @ beta_feature
        resid_outcome_train = train[OUTCOME_COLUMN].to_numpy(dtype=float) - x_train_conf @ beta_outcome
        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            resid_pred = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            resid_feature_test = test[feature_col].to_numpy(dtype=float) - x_test_conf @ beta_feature
            resid_pred = float(slope * resid_feature_test[0])

        resid_actual = float(
            test[OUTCOME_COLUMN].to_numpy(dtype=float)[0] - (x_test_conf @ beta_outcome)[0]
        )
        raw_pred = float((x_test_conf @ beta_outcome)[0] + resid_pred) if math.isfinite(resid_pred) else float("nan")
        raw_outcome = float(test[OUTCOME_COLUMN].iloc[0])
        feature_value = float(test[feature_col].iloc[0])

        rows.append(
            {
                "donor_id": str(test["donor_id"].iloc[0]),
                "outcome": raw_outcome,
                "predicted": raw_pred,
                "actual_residualized_outcome": resid_actual,
                "predicted_residualized_outcome": resid_pred,
                feature_col: feature_value,
            }
        )
        resid_actuals.append(resid_actual)
        resid_preds.append(resid_pred)

    loo_r = _safe_corr(np.asarray(resid_preds, dtype=float), np.asarray(resid_actuals, dtype=float))
    return loo_r, pd.DataFrame(rows)


def _compute_adjusted_metrics(partial_r: float, n_analyzable: int, n_total: int, loo_predictive_r: float) -> dict[str, float]:
    coverage = float(n_analyzable) / float(n_total) if n_total else float("nan")
    selection_score = abs(float(partial_r)) * coverage if math.isfinite(float(partial_r)) else float("nan")
    abs_loo = abs(float(loo_predictive_r)) if math.isfinite(float(loo_predictive_r)) else float("nan")
    is_loo_gap = selection_score - abs_loo if math.isfinite(selection_score) and math.isfinite(abs_loo) else float("nan")
    penalty = max(0.0, is_loo_gap) if math.isfinite(is_loo_gap) else float("nan")
    adjusted_score = selection_score - penalty if math.isfinite(selection_score) and math.isfinite(penalty) else float("nan")
    return {
        "coverage_ratio": coverage,
        "selection_score": selection_score,
        "is_loo_gap": is_loo_gap,
        "penalty": penalty,
        "adjusted_score": adjusted_score,
    }


def _query_any_neighbor(tree: cKDTree | None, points: np.ndarray, radius: float) -> np.ndarray:
    if points.shape[0] == 0:
        return np.zeros(0, dtype=bool)
    if tree is None:
        return np.zeros(points.shape[0], dtype=bool)
    neighbor_lists = tree.query_ball_point(points, r=radius, return_sorted=False)
    return np.fromiter((len(lst) > 0 for lst in neighbor_lists), count=points.shape[0], dtype=bool)


def extract_donor_family_features(
    *,
    donor_id: str,
    data_root: str | Path,
) -> DonorFeatureBundle:
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return DonorFeatureBundle(features={v["feature_column"]: float("nan") for v in VARIATIONS}, counts={})
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name

    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    cells = cells.loc[:, ["x", "y", "cell_type", "region"]].copy()
    cells = cells.loc[cells["region"] == "CA1"].reset_index(drop=True)

    pyramidal = cells.loc[cells["cell_type"] == "Pyramidal Neuron", ["x", "y"]].to_numpy(dtype=float)
    reactive = cells.loc[cells["cell_type"] == "Reactive Astrocyte", ["x", "y"]].to_numpy(dtype=float)
    lymphocyte = cells.loc[cells["cell_type"] == "Lymphocyte", ["x", "y"]].to_numpy(dtype=float)

    counts = {
        "ca1_total_classified_count": int(len(cells)),
        "ca1_pyramidal_count": int(len(pyramidal)),
        "ca1_reactive_astrocyte_count": int(len(reactive)),
        "ca1_lymphocyte_count": int(len(lymphocyte)),
    }

    if len(pyramidal) == 0:
        features = {v["feature_column"]: float("nan") for v in VARIATIONS}
        return DonorFeatureBundle(features=features, counts=counts)

    reactive_tree = cKDTree(reactive) if len(reactive) else None
    lymph_tree = cKDTree(lymphocyte) if len(lymphocyte) else None

    features: dict[str, float] = {}
    for variation in VARIATIONS:
        radius = float(variation["radius_px"])
        has_reactive = _query_any_neighbor(reactive_tree, pyramidal, radius)
        has_lymph = _query_any_neighbor(lymph_tree, pyramidal, radius)
        qualifies = has_reactive & has_lymph
        features[variation["feature_column"]] = float(np.mean(qualifies)) if len(qualifies) else float("nan")

    return DonorFeatureBundle(features=features, counts=counts)


def _load_canonical_variation_name() -> str:
    if BEST_VARIATION_SIDECAR.exists():
        try:
            payload = json.loads(BEST_VARIATION_SIDECAR.read_text())
            name = payload.get("best_variation")
            if name in VARIATION_BY_NAME:
                return str(name)
        except Exception:
            pass
    return DEFAULT_CANONICAL_VARIATION


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Return the winning round-10 biomarker score for one donor from raw data.
    """
    variation_name = _load_canonical_variation_name()
    bundle = extract_donor_family_features(donor_id=donor_id, data_root=data_root)
    feature_column = VARIATION_BY_NAME[variation_name]["feature_column"]
    value = bundle.features.get(feature_column, float("nan"))
    return None if not math.isfinite(float(value)) else float(value)


def main() -> None:
    data_root = Path("/data")
    scratch = Path("/scratch")
    cohort = load_training_cohort(data_root).copy()

    cohort["sex_binary"] = cohort["sex"].map({"Female": 0.0, "Male": 1.0}).astype(float)

    donor_rows: list[dict[str, object]] = []
    for _, row in cohort.iterrows():
        donor_id = str(row["donor_id"])
        bundle = extract_donor_family_features(donor_id=donor_id, data_root=data_root)
        record = row.to_dict()
        record.update(bundle.features)
        record.update(bundle.counts)
        donor_rows.append(record)

    donor_df = pd.DataFrame(donor_rows)

    ranked_variations: list[dict[str, object]] = []
    loo_tables: dict[str, pd.DataFrame] = {}

    for variation in VARIATIONS:
        feature_col = variation["feature_column"]
        full = partial_correlation(
            donor_df,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUNDS,
        )
        analyzable = _clean_analysis_frame(donor_df, feature_col=feature_col)
        n_analyzable = int(len(analyzable))
        loo_r, loo_table = _residualized_loo_predictions(donor_df, feature_col=feature_col)
        loo_tables[variation["variation_name"]] = loo_table
        loo_diag = leave_one_out_summary(
            donor_df,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUNDS,
            id_col="donor_id",
        )

        metrics = _compute_adjusted_metrics(
            partial_r=_score_to_float(full["partial_r"]),
            n_analyzable=n_analyzable,
            n_total=len(donor_df),
            loo_predictive_r=loo_r,
        )
        ranked_variations.append(
            {
                "variation_name": variation["variation_name"],
                "description": variation["description"],
                "feature_column": feature_col,
                "radius_px": variation["radius_px"],
                "n_total": int(len(donor_df)),
                "n_analyzable": n_analyzable,
                "partial_r": _score_to_float(full["partial_r"]),
                "p_value": _score_to_float(full["p_value"]),
                "selection_score": metrics["selection_score"],
                "loo_predictive_r": loo_r,
                "is_loo_gap": metrics["is_loo_gap"],
                "penalty": metrics["penalty"],
                "adjusted_score": metrics["adjusted_score"],
                "coverage_ratio": metrics["coverage_ratio"],
                "loo_unstable_count": int(loo_diag.get("unstable_count", 0) or 0),
                "loo_max_shift": _score_to_float(loo_diag.get("max_shift")),
                "donor_ids_used": analyzable["donor_id"].astype(str).tolist(),
            }
        )

    ranked_variations.sort(
        key=lambda d: (
            -(-np.inf if not math.isfinite(float(d["selection_score"])) else float(d["selection_score"])),
            -(-np.inf if not math.isfinite(float(d["adjusted_score"])) else float(d["adjusted_score"])),
            -(-np.inf if not math.isfinite(float(d["loo_predictive_r"])) else abs(float(d["loo_predictive_r"]))),
        )
    )
    best = ranked_variations[0]
    best_variation = str(best["variation_name"])
    best_feature_column = str(best["feature_column"])

    BEST_VARIATION_SIDECAR.write_text(
        json.dumps(
            {
                "hypothesis_family": HYPOTHESIS_FAMILY,
                "best_variation": best_variation,
                "feature_column": best_feature_column,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    donor_feature_table_path = scratch / "donor_feature_table.csv"
    write_donor_feature_table(
        donor_feature_table_path,
        donor_df,
        feature_column=best_feature_column,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUNDS,
        id_columns=["donor_id", "slide_name"],
        extra_columns=[
            v["feature_column"] for v in VARIATIONS if v["feature_column"] != best_feature_column
        ]
        + [
            "ca1_total_classified_count",
            "ca1_pyramidal_count",
            "ca1_reactive_astrocyte_count",
            "ca1_lymphocyte_count",
        ],
    )

    results_payload = build_results_payload(
        status="success",
        feature_name=HYPOTHESIS_FAMILY,
        outcome=OUTCOME_COLUMN,
        n_total=int(best["n_total"]),
        n_analyzable=int(best["n_analyzable"]),
        partial_r=float(best["partial_r"]),
        ci_lo=float("nan"),
        ci_hi=float("nan"),
        p_value=float(best["p_value"]),
        loo_predictive_r=float(best["loo_predictive_r"]),
        loo_unstable_count=int(best["loo_unstable_count"]),
        loo_max_shift=float(best["loo_max_shift"]),
        donor_ids_used=list(best["donor_ids_used"]),
        covariates=CONFOUNDS,
        recomputed_from_raw=True,
        registry_written=False,
        artifacts={"donor_feature_table": str(donor_feature_table_path)},
        best_variation=best_variation,
        ranked_variations=ranked_variations,
        feature_column=best_feature_column,
        selection_score=float(best["selection_score"]),
        is_loo_gap=float(best["is_loo_gap"]),
        penalty=float(best["penalty"]),
        adjusted_score=float(best["adjusted_score"]),
    )
    write_results_payload(scratch / "results.json", results_payload)

    best_loo_table = loo_tables[best_variation].copy()
    merge_cols = [
        "donor_id",
        best_feature_column,
        "ca1_pyramidal_count",
        "ca1_reactive_astrocyte_count",
        "ca1_lymphocyte_count",
    ]
    best_loo_table = best_loo_table.merge(
        donor_df.loc[:, merge_cols],
        on="donor_id",
        how="left",
        suffixes=("", "_dup"),
    )

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best_variation}")
    print(f"  IS partial r:      {float(best['partial_r']):.4f}")
    print(f"  Selection score:   {float(best['selection_score']):.4f}")
    print(f"  LOO predictive r:  {float(best['loo_predictive_r']):.4f}  (diagnostic)")
    print(f"  IS-LOO Gap:        {float(best['is_loo_gap']):.4f}  (penalty={float(best['penalty']):.4f})")
    print(f"  Adjusted Score:    {float(best['adjusted_score']):.4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked_variations:
        print(
            f"  {item['variation_name']}  "
            f"{float(item['partial_r']):.4f}  "
            f"{float(item['selection_score']):.4f}  "
            f"{float(item['loo_predictive_r']):.4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(
        f"  donor_id  outcome  predicted  {best_feature_column}  "
        "ca1_pyramidal_count  ca1_reactive_astrocyte_count  ca1_lymphocyte_count"
    )
    for _, row in best_loo_table.sort_values("donor_id").iterrows():
        print(
            f"  {row['donor_id']}  "
            f"{float(row['outcome']):.4f}  "
            f"{float(row['predicted']):.4f}  "
            f"{float(row[best_feature_column]):.6f}  "
            f"{int(row['ca1_pyramidal_count'])}  "
            f"{int(row['ca1_reactive_astrocyte_count'])}  "
            f"{int(row['ca1_lymphocyte_count'])}"
        )


if __name__ == "__main__":
    main()
