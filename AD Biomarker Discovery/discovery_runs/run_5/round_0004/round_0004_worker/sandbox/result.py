from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path
from typing import Any

sys.path.append("/shared/lib")

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort
from shared_analysis.stats import (
    bootstrap_partial_correlation,
    bootstrap_partial_correlation_stability,
    leave_one_out_summary,
    partial_correlation,
    residualized_loo_predictive_correlation,
)

try:
    import openslide  # type: ignore
except Exception:  # pragma: no cover
    openslide = None


HYPOTHESIS_FAMILY = "CA1 peripyramidal reactive astrocyte hypertrophy area"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
DISTANCE_THRESHOLD_UM = 30.0
DEFAULT_MPP_UM = 0.503

VARIATIONS = {
    "candidate_variant_a": {
        "description": (
            "Median area of CA1 reactive astrocytes whose centroids are within "
            "30 microns of any CA1 pyramidal neuron."
        ),
        "feature_column": "ca1_peripyramidal_reactive_astro_area_median_30um_um2",
        "summary": "median",
    },
    "candidate_variant_b": {
        "description": (
            "75th percentile area of CA1 reactive astrocytes whose centroids are within "
            "30 microns of any CA1 pyramidal neuron."
        ),
        "feature_column": "ca1_peripyramidal_reactive_astro_area_p75_30um_um2",
        "summary": "p75",
    },
}

FEATURE_NAME = "ca1_peripyramidal_reactive_astro_area_median_30um_um2"
FEATURE_COLUMN = "ca1_peripyramidal_reactive_astro_area_median_30um_um2"
BEST_VARIATION_FALLBACK = "candidate_variant_a"
BEST_SPEC_PATH = Path(__file__).with_name("best_variation.json")
DONOR_TABLE_PATH = Path(__file__).with_name("donor_feature_table.csv")
RESULTS_PATH = Path(__file__).with_name("results.json")


def _safe_float(value: Any) -> float:
    try:
        out = float(value)
    except Exception:
        return float("nan")
    return out


def _is_finite(value: Any) -> bool:
    try:
        return bool(np.isfinite(float(value)))
    except Exception:
        return False


def _load_best_spec() -> dict[str, Any]:
    if BEST_SPEC_PATH.exists():
        try:
            return json.loads(BEST_SPEC_PATH.read_text(encoding="utf-8"))
        except Exception:
            pass
    fallback = VARIATIONS[BEST_VARIATION_FALLBACK]
    return {
        "best_variation": BEST_VARIATION_FALLBACK,
        "feature_name": fallback["feature_column"],
        "feature_column": fallback["feature_column"],
    }


def _resolve_best_feature_column() -> str:
    spec = _load_best_spec()
    return str(spec.get("feature_column") or FEATURE_COLUMN)


def _resolve_best_feature_name() -> str:
    spec = _load_best_spec()
    return str(spec.get("feature_name") or FEATURE_NAME)


def _sex_to_binary(series: pd.Series) -> pd.Series:
    mapped = series.astype(str).str.strip().str.lower().map({"female": 0.0, "male": 1.0})
    return mapped.astype(float)


def _design_matrix(confounds: pd.DataFrame) -> np.ndarray:
    x = confounds.astype(float).to_numpy()
    intercept = np.ones((len(confounds), 1), dtype=float)
    return np.hstack([intercept, x])


def _clean_frame(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    frame = df.loc[:, columns].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    return frame


def _get_mpp_um(svs_path: Path) -> tuple[float, float]:
    if openslide is not None:
        try:
            slide = openslide.OpenSlide(str(svs_path))
            px = slide.properties.get("openslide.mpp-x") or slide.properties.get("aperio.MPP")
            py = slide.properties.get("openslide.mpp-y") or slide.properties.get("aperio.MPP")
            mpp_x = float(px)
            mpp_y = float(py)
            if np.isfinite(mpp_x) and np.isfinite(mpp_y) and mpp_x > 0 and mpp_y > 0:
                return mpp_x, mpp_y
        except Exception:
            pass
    return DEFAULT_MPP_UM, DEFAULT_MPP_UM


def _extract_single_slide_features(zarr_path: Path, svs_path: Path) -> dict[str, Any]:
    warnings.filterwarnings("ignore", message="Object at .* not recognized as a component of a Zarr hierarchy.")
    warnings.filterwarnings("ignore", message="invalid value encountered in divide")
    cells = build_cell_table(zarr_path, include_regions=True, include_geometry=True)

    region = cells["region"].astype(object)
    cell_type = cells["cell_type"].astype(str)

    ca1_mask = region == "CA1"
    reactive_mask = ca1_mask & (cell_type == "Reactive Astrocyte")
    pyramidal_mask = ca1_mask & (cell_type == "Pyramidal Neuron")

    reactive = cells.loc[reactive_mask, ["x", "y", "area"]].copy()
    pyramidal = cells.loc[pyramidal_mask, ["x", "y"]].copy()

    mpp_x, mpp_y = _get_mpp_um(svs_path)
    distance_threshold_px = DISTANCE_THRESHOLD_UM / float((mpp_x + mpp_y) / 2.0)
    area_scale = float(mpp_x * mpp_y)

    qualifying_areas_um2: np.ndarray
    if len(reactive) == 0 or len(pyramidal) == 0:
        qualifying_areas_um2 = np.asarray([], dtype=float)
    else:
        tree = cKDTree(pyramidal[["x", "y"]].to_numpy(dtype=float))
        dists, _ = tree.query(
            reactive[["x", "y"]].to_numpy(dtype=float),
            k=1,
            workers=-1,
        )
        dists = np.asarray(dists, dtype=float)
        near_mask = dists <= distance_threshold_px
        qualifying_areas_um2 = reactive["area"].to_numpy(dtype=float)[near_mask] * area_scale

    median_area = float(np.nanmedian(qualifying_areas_um2)) if len(qualifying_areas_um2) else float("nan")
    p75_area = float(np.nanpercentile(qualifying_areas_um2, 75)) if len(qualifying_areas_um2) else float("nan")

    return {
        VARIATIONS["candidate_variant_a"]["feature_column"]: median_area,
        VARIATIONS["candidate_variant_b"]["feature_column"]: p75_area,
        "ca1_reactive_astro_count": int(len(reactive)),
        "ca1_pyramidal_neuron_count": int(len(pyramidal)),
        "ca1_peripyramidal_reactive_astro_count_30um": int(len(qualifying_areas_um2)),
        "ca1_peripyramidal_reactive_astro_area_mean_30um_um2": (
            float(np.mean(qualifying_areas_um2)) if len(qualifying_areas_um2) else float("nan")
        ),
        "ca1_peripyramidal_reactive_astro_area_iqr_30um_um2": (
            float(np.subtract(*np.percentile(qualifying_areas_um2, [75, 25])))
            if len(qualifying_areas_um2) >= 2
            else float("nan")
        ),
        "distance_threshold_um": float(DISTANCE_THRESHOLD_UM),
        "distance_threshold_px": float(distance_threshold_px),
        "mpp_x_um": float(mpp_x),
        "mpp_y_um": float(mpp_y),
    }


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Canonical replay target: returns the winning variation saved in best_variation.json.
    """
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    zarr_path = data_root / slide_name
    svs_path = zarr_path.with_suffix("")
    features = _extract_single_slide_features(zarr_path, svs_path)
    feature_column = _resolve_best_feature_column()
    return _safe_float(features.get(feature_column, np.nan))


def _build_feature_table(data_root: Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = _sex_to_binary(cohort["sex"])

    rows: list[dict[str, Any]] = []
    for _, row in cohort.iterrows():
        donor_id = str(row["donor_id"])
        slide_name = str(row["slide_name"])
        zarr_path = data_root / slide_name
        svs_path = zarr_path.with_suffix("")
        feature_row = {
            "donor_id": donor_id,
            "slide_name": slide_name,
            OUTCOME_COLUMN: _safe_float(row[OUTCOME_COLUMN]),
            "max_age_vis": _safe_float(row["max_age_vis"]),
            "braak_numeric": _safe_float(row["braak_numeric"]),
            "cerad_ordinal": _safe_float(row["cerad_ordinal"]),
            "sex": row["sex"],
            "sex_binary": _safe_float(row["sex_binary"]),
        }
        feature_row.update(_extract_single_slide_features(zarr_path, svs_path))
        rows.append(feature_row)

    table = pd.DataFrame(rows)
    ordered_cols = [
        "donor_id",
        "slide_name",
        OUTCOME_COLUMN,
        "max_age_vis",
        "braak_numeric",
        "cerad_ordinal",
        "sex",
        "sex_binary",
        VARIATIONS["candidate_variant_a"]["feature_column"],
        VARIATIONS["candidate_variant_b"]["feature_column"],
        "ca1_reactive_astro_count",
        "ca1_pyramidal_neuron_count",
        "ca1_peripyramidal_reactive_astro_count_30um",
        "ca1_peripyramidal_reactive_astro_area_mean_30um_um2",
        "ca1_peripyramidal_reactive_astro_area_iqr_30um_um2",
        "distance_threshold_um",
        "distance_threshold_px",
        "mpp_x_um",
        "mpp_y_um",
    ]
    table = table.loc[:, ordered_cols]
    table.to_csv(DONOR_TABLE_PATH, index=False)
    return table


def _loo_prediction_table(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
    donor_col: str = "donor_id",
) -> pd.DataFrame:
    cols = [donor_col, feature_col, outcome_col, *confounds]
    frame = _clean_frame(df, cols)
    rows: list[dict[str, Any]] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train_conf = _design_matrix(train[confounds])
        x_test_conf = _design_matrix(test[confounds])

        y_feat_train = train[feature_col].to_numpy(dtype=float)
        y_out_train = train[outcome_col].to_numpy(dtype=float)

        beta_feature, *_ = np.linalg.lstsq(x_train_conf, y_feat_train, rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train_conf, y_out_train, rcond=None)

        resid_feature_train = y_feat_train - x_train_conf @ beta_feature
        resid_outcome_train = y_out_train - x_train_conf @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0 or not np.isfinite(denom):
            pred_resid = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            test_feat = test[feature_col].to_numpy(dtype=float)
            pred_feat_test = x_test_conf @ beta_feature
            pred_resid = float(slope * (test_feat[0] - pred_feat_test[0]))

        actual_outcome = float(test[outcome_col].iloc[0])
        base_outcome = float((x_test_conf @ beta_outcome)[0])
        predicted_outcome = base_outcome + pred_resid if np.isfinite(pred_resid) else float("nan")
        actual_resid = actual_outcome - base_outcome
        abs_error = abs(predicted_outcome - actual_outcome) if np.isfinite(predicted_outcome) else float("nan")

        rows.append(
            {
                donor_col: str(test[donor_col].iloc[0]),
                "outcome": actual_outcome,
                "predicted": predicted_outcome,
                "predicted_residual": pred_resid,
                "actual_residual": actual_resid,
                feature_col: float(test[feature_col].iloc[0]),
                "abs_error": abs_error,
            }
        )

    return pd.DataFrame(rows)


def _evaluate_variation(feature_table: pd.DataFrame, variation_name: str) -> dict[str, Any]:
    spec = VARIATIONS[variation_name]
    feature_col = str(spec["feature_column"])

    base_metrics = partial_correlation(
        feature_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    seed_base = {"candidate_variant_a": 101, "candidate_variant_b": 202}[variation_name]
    boot_metrics = bootstrap_partial_correlation(
        feature_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        n_boot=2000,
        random_state=seed_base,
    )
    stability_metrics = bootstrap_partial_correlation_stability(
        feature_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        n_boot=400,
        sample_frac=0.8,
        random_state=seed_base + 17,
    )
    loo_predictive_r = residualized_loo_predictive_correlation(
        feature_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    loo_summary = leave_one_out_summary(
        feature_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        id_col="donor_id",
        unstable_delta=0.10,
    )

    n_total = int(len(feature_table))
    n_analyzable = int(base_metrics["n"])
    coverage_ratio = float(n_analyzable / n_total) if n_total else float("nan")
    partial_r_value = _safe_float(base_metrics["partial_r"])
    selection_score = abs(partial_r_value) * coverage_ratio if _is_finite(partial_r_value) else float("nan")
    gap = (
        abs(partial_r_value - float(loo_predictive_r))
        if _is_finite(partial_r_value) and _is_finite(loo_predictive_r)
        else float("nan")
    )
    penalty = gap if _is_finite(gap) else float("nan")
    adjusted_score = (
        selection_score - penalty if _is_finite(selection_score) and _is_finite(penalty) else float("nan")
    )
    per_donor_loo = _loo_prediction_table(
        feature_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        donor_col="donor_id",
    )

    donor_ids_used = (
        _clean_frame(feature_table, ["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS])["donor_id"]
        .astype(str)
        .tolist()
    )

    result = {
        "name": variation_name,
        "description": spec["description"],
        "feature_column": feature_col,
        "feature_name": feature_col,
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "coverage_ratio": coverage_ratio,
        "partial_r": partial_r_value,
        "p_value": _safe_float(base_metrics["p_value"]),
        "ci_lo": _safe_float(boot_metrics.get("ci_lo", np.nan)),
        "ci_hi": _safe_float(boot_metrics.get("ci_hi", np.nan)),
        "selection_score": selection_score,
        "loo_predictive_r": _safe_float(loo_predictive_r),
        "gap": gap,
        "penalty": penalty,
        "adjusted_score": adjusted_score,
        "loo_unstable_count": int(loo_summary.get("unstable_count", 0) or 0),
        "loo_max_shift": _safe_float(loo_summary.get("max_shift", np.nan)),
        "unstable_donors": loo_summary.get("unstable_donors", []),
        "bootstrap_median_partial_r": _safe_float(stability_metrics.get("bootstrap_median_partial_r", np.nan)),
        "bootstrap_q25_partial_r": _safe_float(stability_metrics.get("bootstrap_q25_partial_r", np.nan)),
        "bootstrap_q75_partial_r": _safe_float(stability_metrics.get("bootstrap_q75_partial_r", np.nan)),
        "bootstrap_sign_consistency": _safe_float(stability_metrics.get("bootstrap_sign_consistency", np.nan)),
        "bootstrap_positive_fraction": _safe_float(stability_metrics.get("bootstrap_positive_fraction", np.nan)),
        "bootstrap_negative_fraction": _safe_float(stability_metrics.get("bootstrap_negative_fraction", np.nan)),
        "donor_ids_used": donor_ids_used,
        "per_donor_loo": per_donor_loo.to_dict(orient="records"),
    }
    return result


def _rank_variations(feature_table: pd.DataFrame) -> list[dict[str, Any]]:
    ranked = [_evaluate_variation(feature_table, name) for name in VARIATIONS]
    ranked.sort(
        key=lambda r: (
            -999.0 if not _is_finite(r["selection_score"]) else -float(r["selection_score"]),
            -999.0 if not _is_finite(r["partial_r"]) else -abs(float(r["partial_r"])),
            -999.0 if not _is_finite(r["loo_predictive_r"]) else -float(r["loo_predictive_r"]),
        )
    )
    return ranked


def _format_num(value: Any) -> str:
    try:
        value = float(value)
    except Exception:
        return "nan"
    if not np.isfinite(value):
        return "nan"
    return f"{value:.4f}"


def _write_best_sidecar(best: dict[str, Any]) -> None:
    payload = {
        "best_variation": best["name"],
        "feature_name": best["feature_name"],
        "feature_column": best["feature_column"],
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "distance_threshold_um": DISTANCE_THRESHOLD_UM,
    }
    BEST_SPEC_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _write_results_json(feature_table: pd.DataFrame, ranked: list[dict[str, Any]]) -> dict[str, Any]:
    best = ranked[0]
    payload = {
        "status": "ok",
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "feature_name": best["feature_name"],
        "feature_column": best["feature_column"],
        "best_variation": best["name"],
        "outcome": OUTCOME_COLUMN,
        "covariates": CONFOUND_COLUMNS,
        "n_total": int(best["n_total"]),
        "n_analyzable": int(best["n_analyzable"]),
        "coverage_ratio": float(best["coverage_ratio"]),
        "partial_r": float(best["partial_r"]),
        "p_value": float(best["p_value"]),
        "ci_lo": float(best["ci_lo"]),
        "ci_hi": float(best["ci_hi"]),
        "selection_score": float(best["selection_score"]),
        "loo_predictive_r": float(best["loo_predictive_r"]),
        "is_loo_gap": float(best["gap"]),
        "gap_penalty": float(best["penalty"]),
        "adjusted_score": float(best["adjusted_score"]),
        "loo_unstable_count": int(best["loo_unstable_count"]),
        "loo_max_shift": float(best["loo_max_shift"]),
        "bootstrap_median_partial_r": float(best["bootstrap_median_partial_r"]),
        "bootstrap_sign_consistency": float(best["bootstrap_sign_consistency"]),
        "ranked_variations": ranked,
        "donor_ids_used": best["donor_ids_used"],
        "recomputed_from_raw": True,
        "registry_written": False,
        "artifacts": {
            "donor_feature_table": str(DONOR_TABLE_PATH),
            "best_variation_sidecar": str(BEST_SPEC_PATH),
        },
        "current_panel_score": 0.5984,
        "accepted_panel_members": [
            "ca1_reactive_astro_fraction",
            "ca1_pyramidal_near_reactive_astro_fraction_30um",
            "ca1_pyramidal_reactive_astro_crowding_fraction_k3_30um",
        ],
    }
    RESULTS_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


def _print_stdout(best: dict[str, Any], ranked: list[dict[str, Any]]) -> None:
    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best['name']}")
    print(f"  IS partial r:      {_format_num(best['partial_r'])}")
    print(f"  Selection score:   {_format_num(best['selection_score'])}")
    print(f"  LOO predictive r:  {_format_num(best['loo_predictive_r'])}  (diagnostic)")
    print(f"  IS-LOO Gap:        {_format_num(best['gap'])}  (penalty={_format_num(best['penalty'])})")
    print(f"  Adjusted Score:    {_format_num(best['adjusted_score'])}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for entry in ranked:
        print(
            "  "
            f"{entry['name']}  "
            f"{_format_num(entry['partial_r'])}  "
            f"{_format_num(entry['selection_score'])}  "
            f"{_format_num(entry['loo_predictive_r'])}"
        )
    print()
    print("PER-DONOR (LOO):")
    feature_col = best["feature_column"]
    print(f"  donor_id  outcome  predicted  {feature_col}")
    for row in best["per_donor_loo"]:
        print(
            "  "
            f"{row['donor_id']}  "
            f"{_format_num(row['outcome'])}  "
            f"{_format_num(row['predicted'])}  "
            f"{_format_num(row[feature_col])}"
        )


def main() -> int:
    data_root = Path("/data")
    feature_table = _build_feature_table(data_root)
    ranked = _rank_variations(feature_table)
    best = ranked[0]
    _write_best_sidecar(best)
    _write_results_json(feature_table, ranked)
    _print_stdout(best, ranked)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
