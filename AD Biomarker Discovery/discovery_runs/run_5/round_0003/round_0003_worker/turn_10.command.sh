cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.append("/shared/lib")

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from shared_analysis import build_cell_table, load_training_cohort
from shared_analysis.stats import partial_correlation

FAMILY_NAME = "ca1_peripyramidal_reactive_astro_crowding"
OUTCOME_COL = "slope_zmem0"
CONFOUNDS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]

RADIUS_UM = 30.0
DEFAULT_MPP_UM_PER_PX = 0.5

VARIATIONS = [
    {
        "name": "k2_within_30um",
        "description": "Fraction of CA1 pyramidal neurons with at least 2 CA1 reactive astrocytes within 30 um.",
        "k_threshold": 2,
        "radius_um": RADIUS_UM,
        "feature_column": "ca1_pyramidal_reactive_astro_crowding_fraction_k2_30um",
    },
    {
        "name": "k3_within_30um",
        "description": "Fraction of CA1 pyramidal neurons with at least 3 CA1 reactive astrocytes within 30 um.",
        "k_threshold": 3,
        "radius_um": RADIUS_UM,
        "feature_column": "ca1_pyramidal_reactive_astro_crowding_fraction_k3_30um",
    },
]

# Canonical replay target for evaluator import. Updated to the winning variation
# after the worker run if needed.
CANONICAL_VARIATION_NAME = "k2_within_30um"
FEATURE_NAME = "ca1_pyramidal_reactive_astro_crowding_fraction_k2_30um"
FEATURE_COLUMN = "ca1_pyramidal_reactive_astro_crowding_fraction_k2_30um"


def _variation_lookup() -> dict[str, dict]:
    return {v["name"]: v for v in VARIATIONS}


def _safe_float(value) -> float | None:
    if value is None:
        return None
    try:
        value = float(value)
    except Exception:
        return None
    if not math.isfinite(value):
        return None
    return value


def _sex_to_binary(series: pd.Series) -> pd.Series:
    return series.map({"Female": 1.0, "Male": 0.0}).astype(float)


def _design_matrix(df: pd.DataFrame) -> np.ndarray:
    x = np.column_stack([np.ones(len(df), dtype=float)] + [df[col].to_numpy(dtype=float) for col in CONFOUNDS])
    return x


def _residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ beta


def _read_slide_mpp_um_per_px(slide_path: Path) -> float:
    svs_path = slide_path.with_suffix("")
    if svs_path.suffix != ".svs":
        svs_path = slide_path.parent / slide_path.name.replace(".zarr", "")
    try:
        import openslide  # type: ignore

        with openslide.OpenSlide(str(svs_path)) as slide:
            props = slide.properties
            mpp_x = props.get(openslide.PROPERTY_NAME_MPP_X) or props.get("aperio.MPP")
            mpp_y = props.get(openslide.PROPERTY_NAME_MPP_Y) or props.get("aperio.MPP")
            vals = [float(v) for v in [mpp_x, mpp_y] if v is not None]
            if vals:
                return float(np.mean(vals))
    except Exception:
        pass
    return DEFAULT_MPP_UM_PER_PX


def _count_neighbors_within_radius(anchor_xy: np.ndarray, neighbor_xy: np.ndarray, radius_px: float) -> np.ndarray:
    if len(anchor_xy) == 0:
        return np.zeros(0, dtype=int)
    if len(neighbor_xy) == 0:
        return np.zeros(len(anchor_xy), dtype=int)
    tree = cKDTree(neighbor_xy)
    try:
        counts = tree.query_ball_point(anchor_xy, r=radius_px, return_length=True)
        return np.asarray(counts, dtype=int)
    except TypeError:
        neighborhoods = tree.query_ball_point(anchor_xy, r=radius_px)
        return np.fromiter((len(v) for v in neighborhoods), dtype=int, count=len(anchor_xy))


def _compute_variation_scores_for_slide(slide_path: Path) -> dict[str, float | None]:
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    ca1 = cells.loc[cells["region"] == "CA1", ["x", "y", "cell_type"]].copy()

    pyramidal = ca1.loc[ca1["cell_type"] == "Pyramidal Neuron", ["x", "y"]].to_numpy(dtype=float)
    reactive = ca1.loc[ca1["cell_type"] == "Reactive Astrocyte", ["x", "y"]].to_numpy(dtype=float)

    if len(pyramidal) == 0:
        return {v["name"]: np.nan for v in VARIATIONS}

    mpp = _read_slide_mpp_um_per_px(slide_path)
    radius_px = float(RADIUS_UM / mpp)
    neighbor_counts = _count_neighbors_within_radius(pyramidal, reactive, radius_px)

    scores: dict[str, float | None] = {}
    for variation in VARIATIONS:
        frac = float(np.mean(neighbor_counts >= int(variation["k_threshold"])))
        scores[variation["name"]] = frac
    return scores


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Return the canonical winning biomarker score for one donor.

    The evaluator imports this function later. It recomputes the same raw-cell
    donor scalar from the .zarr slide and uses the canonical winning variation.
    """
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    scores = _compute_variation_scores_for_slide(data_root / slide_name)
    value = scores.get(CANONICAL_VARIATION_NAME)
    return None if value is None or not math.isfinite(float(value)) else float(value)


def _residualized_loo_table(df: pd.DataFrame, feature_col: str) -> tuple[float, pd.DataFrame]:
    frame = df.dropna(subset=[feature_col, OUTCOME_COL, *CONFOUNDS]).reset_index(drop=True).copy()
    rows: list[dict] = []
    preds: list[float] = []
    actuals: list[float] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train_conf = _design_matrix(train)
        x_test_conf = _design_matrix(test)

        feat_train = train[feature_col].to_numpy(dtype=float)
        y_train = train[OUTCOME_COL].to_numpy(dtype=float)

        beta_feature, *_ = np.linalg.lstsq(x_train_conf, feat_train, rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train_conf, y_train, rcond=None)

        resid_feature_train = feat_train - x_train_conf @ beta_feature
        resid_outcome_train = y_train - x_train_conf @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            pred_resid = float("nan")
            actual_resid = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            feat_test = test[feature_col].to_numpy(dtype=float)
            y_test = test[OUTCOME_COL].to_numpy(dtype=float)
            resid_feature_test = feat_test - x_test_conf @ beta_feature
            actual_resid = float(y_test[0] - (x_test_conf @ beta_outcome)[0])
            pred_resid = float(slope * resid_feature_test[0])

        preds.append(pred_resid)
        actuals.append(actual_resid)
        rows.append(
            {
                "donor_id": test.iloc[0]["donor_id"],
                "outcome": float(test.iloc[0][OUTCOME_COL]),
                "actual_residualized": actual_resid,
                "predicted": pred_resid,
                feature_col: float(test.iloc[0][feature_col]),
                "abs_error": abs(actual_resid - pred_resid) if math.isfinite(actual_resid) and math.isfinite(pred_resid) else np.nan,
            }
        )

    pred_arr = np.asarray(preds, dtype=float)
    act_arr = np.asarray(actuals, dtype=float)
    mask = np.isfinite(pred_arr) & np.isfinite(act_arr)
    if mask.sum() < 3 or np.std(pred_arr[mask]) == 0 or np.std(act_arr[mask]) == 0:
        corr = float("nan")
    else:
        corr = float(np.corrcoef(pred_arr[mask], act_arr[mask])[0, 1])
    return corr, pd.DataFrame(rows)


def _compute_metrics(feature_df: pd.DataFrame, variation: dict) -> tuple[dict, pd.DataFrame]:
    feature_col = variation["feature_column"]
    analyzable = feature_df.dropna(subset=[feature_col, OUTCOME_COL, *CONFOUNDS]).copy()
    n_total = int(len(feature_df))
    n_analyzable = int(len(analyzable))
    coverage_ratio = float(n_analyzable / n_total) if n_total else 0.0

    pc = partial_correlation(
        feature_df,
        feature_col=feature_col,
        outcome_col=OUTCOME_COL,
        confounds=CONFOUNDS,
    )
    partial_r = float(pc["partial_r"]) if pc.get("partial_r") is not None else float("nan")
    p_value = float(pc["p_value"]) if pc.get("p_value") is not None else float("nan")
    selection_score = float(abs(partial_r) * coverage_ratio) if math.isfinite(partial_r) else float("nan")

    loo_r, loo_table = _residualized_loo_table(feature_df, feature_col=feature_col)
    is_loo_gap = float(abs(partial_r) - loo_r) if math.isfinite(partial_r) and math.isfinite(loo_r) else float("nan")
    gap_penalty = float(max(is_loo_gap, 0.0)) if math.isfinite(is_loo_gap) else float("nan")
    adjusted_score = float(selection_score - gap_penalty) if math.isfinite(selection_score) and math.isfinite(gap_penalty) else float("nan")

    metrics = {
        "variation_name": variation["name"],
        "description": variation["description"],
        "radius_um": variation["radius_um"],
        "k_threshold": variation["k_threshold"],
        "feature_name": feature_col,
        "feature_column": feature_col,
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "coverage_ratio": coverage_ratio,
        "partial_r": partial_r,
        "p_value": p_value,
        "selection_score": selection_score,
        "loo_predictive_r": loo_r,
        "is_loo_gap": is_loo_gap,
        "gap_penalty": gap_penalty,
        "adjusted_score": adjusted_score,
    }
    return metrics, loo_table


def _format_float(value: float | None) -> str:
    if value is None:
        return "nan"
    try:
        value = float(value)
    except Exception:
        return "nan"
    if not math.isfinite(value):
        return "nan"
    return f"{value:.4f}"


def main() -> int:
    data_root = Path("/data")
    scratch_dir = Path("/scratch")
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = _sex_to_binary(cohort["sex"])

    donor_rows: list[dict] = []
    for _, row in cohort.iterrows():
        slide_path = data_root / str(row["slide_name"])
        scores = _compute_variation_scores_for_slide(slide_path)
        out_row = {
            "donor_id": row["donor_id"],
            "slide_name": row["slide_name"],
            OUTCOME_COL: row[OUTCOME_COL],
            "max_age_vis": row["max_age_vis"],
            "braak_numeric": row["braak_numeric"],
            "cerad_ordinal": row["cerad_ordinal"],
            "sex_binary": row["sex_binary"],
        }
        for variation in VARIATIONS:
            out_row[variation["feature_column"]] = scores.get(variation["name"])
        donor_rows.append(out_row)

    feature_df = pd.DataFrame(donor_rows)
    feature_df.to_csv(scratch_dir / "donor_feature_table.csv", index=False)

    ranked: list[dict] = []
    loo_tables: dict[str, pd.DataFrame] = {}
    for variation in VARIATIONS:
        metrics, loo_table = _compute_metrics(feature_df, variation)
        ranked.append(metrics)
        loo_tables[variation["name"]] = loo_table

    ranked.sort(
        key=lambda d: (
            -999 if not math.isfinite(float(d["selection_score"])) else -float(d["selection_score"]),
            -999 if not math.isfinite(float(d["loo_predictive_r"])) else -float(d["loo_predictive_r"]),
        )
    )
    best = ranked[0]
    best_variation_name = str(best["variation_name"])
    best_feature_col = str(best["feature_column"])
    best_loo_table = loo_tables[best_variation_name].copy()

    results_payload = {
        "hypothesis_family": FAMILY_NAME,
        "feature_name": best_feature_col,
        "feature_column": best_feature_col,
        "best_variation": best_variation_name,
        "partial_r": best["partial_r"],
        "selection_score": best["selection_score"],
        "loo_predictive_r": best["loo_predictive_r"],
        "is_loo_gap": best["is_loo_gap"],
        "gap_penalty": best["gap_penalty"],
        "adjusted_score": best["adjusted_score"],
        "n_total": best["n_total"],
        "n_analyzable": best["n_analyzable"],
        "coverage_ratio": best["coverage_ratio"],
        "p_value": best["p_value"],
        "ranked_variations": ranked,
    }
    (scratch_dir / "results.json").write_text(json.dumps(results_payload, indent=2))

    print(f"HYPOTHESIS FAMILY: {FAMILY_NAME}")
    print(f"BEST VARIATION: {best_variation_name}")
    print(f"  IS partial r:      {_format_float(best['partial_r'])}")
    print(f"  Selection score:   {_format_float(best['selection_score'])}")
    print(f"  LOO predictive r:  {_format_float(best['loo_predictive_r'])}  (diagnostic)")
    print(f"  IS-LOO Gap:        {_format_float(best['is_loo_gap'])}  (penalty={_format_float(best['gap_penalty'])})")
    print(f"  Adjusted Score:    {_format_float(best['adjusted_score'])}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked:
        print(
            f"  {item['variation_name']}  "
            f"{_format_float(item['partial_r'])}  "
            f"{_format_float(item['selection_score'])}  "
            f"{_format_float(item['loo_predictive_r'])}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best_feature_col}")
    for _, row in best_loo_table.sort_values("donor_id").iterrows():
        print(
            f"  {row['donor_id']}  "
            f"{_format_float(row['outcome'])}  "
            f"{_format_float(row['predicted'])}  "
            f"{_format_float(row[best_feature_col])}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
PY
python /scratch/result.py > /scratch/run_stdout.txt
cat /scratch/run_stdout.txt
