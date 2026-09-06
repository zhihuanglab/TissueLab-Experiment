cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path

SHARED_LIB = Path("/shared/lib")
if SHARED_LIB.exists() and str(SHARED_LIB) not in sys.path:
    sys.path.insert(0, str(SHARED_LIB))

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from shared_analysis import build_cell_table, load_training_cohort
from shared_analysis.artifacts import (
    build_results_payload,
    write_donor_feature_table,
    write_results_payload,
)
from shared_analysis.stats import partial_correlation

DATA_ROOT_DEFAULT = Path("/data")
OUTCOME_COLUMN = "slope_zmem0"
CONFOUNDS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
FEATURE_FAMILY = "ca1_reactive_astrocyte_proximal_pyramidal_area"
BASELINE_VARIATION = "median_log_area_60px"
MIN_SELECTED_CELLS = 10

VARIATIONS = {
    "median_log_area_60px": {
        "radius_px": 60.0,
        "feature_column": "ca1_pyramidal_reactive_astro_proximal_median_log_area_60px",
        "selected_count_column": "ca1_pyramidal_reactive_astro_proximal_count_60px",
        "description": (
            "Median log1p(contour area px^2) of CA1 pyramidal neurons within 60 px "
            "of any CA1 reactive astrocyte centroid."
        ),
    },
    "median_log_area_80px": {
        "radius_px": 80.0,
        "feature_column": "ca1_pyramidal_reactive_astro_proximal_median_log_area_80px",
        "selected_count_column": "ca1_pyramidal_reactive_astro_proximal_count_80px",
        "description": (
            "Median log1p(contour area px^2) of CA1 pyramidal neurons within 80 px "
            "of any CA1 reactive astrocyte centroid."
        ),
    },
}
KEY_COUNT_COLUMNS = ["ca1_pyramidal_neuron_count", "ca1_reactive_astrocyte_count"] + [
    spec["selected_count_column"] for spec in VARIATIONS.values()
]


def _safe_float(value):
    if value is None:
        return None
    try:
        value = float(value)
    except Exception:
        return None
    return value if math.isfinite(value) else None


def _fmt(value: float | None) -> str:
    if value is None:
        return "nan"
    try:
        value = float(value)
    except Exception:
        return "nan"
    if not math.isfinite(value):
        return "nan"
    return f"{value:.4f}"


def _results_sidecar_path() -> Path:
    return Path(__file__).with_name("results.json")


def _load_best_variation_from_sidecar() -> str:
    sidecar = _results_sidecar_path()
    if sidecar.exists():
        try:
            payload = json.loads(sidecar.read_text(encoding="utf-8"))
            name = payload.get("best_variation")
            if name in VARIATIONS:
                return str(name)
        except Exception:
            pass
    return BASELINE_VARIATION


CANONICAL_VARIATION = _load_best_variation_from_sidecar()
FEATURE_NAME = VARIATIONS[CANONICAL_VARIATION]["feature_column"]
FEATURE_COLUMN = FEATURE_NAME


def _prepare_cohort(data_root: str | Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = (cohort["sex"].astype(str).str.lower() == "male").astype(float)
    return cohort


def _design_matrix(confounds: pd.DataFrame) -> np.ndarray:
    matrix = confounds.astype(float).to_numpy()
    intercept = np.ones((len(confounds), 1), dtype=float)
    return np.hstack([intercept, matrix])


def _extract_slide_features(slide_path: Path) -> dict[str, float]:
    warnings.filterwarnings("ignore", category=UserWarning)
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=True)

    mask = (cells["region"] == "CA1") & (
        cells["cell_type"].isin(["Pyramidal Neuron", "Reactive Astrocyte"])
    )
    ca1 = cells.loc[mask, ["x", "y", "cell_type", "area"]].copy()

    pyramidal = ca1.loc[ca1["cell_type"] == "Pyramidal Neuron", ["x", "y", "area"]].copy()
    reactive = ca1.loc[ca1["cell_type"] == "Reactive Astrocyte", ["x", "y"]].copy()

    result: dict[str, float] = {
        "ca1_pyramidal_neuron_count": float(len(pyramidal)),
        "ca1_reactive_astrocyte_count": float(len(reactive)),
    }
    for spec in VARIATIONS.values():
        result[spec["selected_count_column"]] = 0.0

    if len(pyramidal) == 0 or len(reactive) == 0:
        for spec in VARIATIONS.values():
            result[spec["feature_column"]] = float("nan")
        return result

    pyr_xy = pyramidal[["x", "y"]].to_numpy(dtype=np.float32, copy=False)
    pyr_log_area = np.log1p(pyramidal["area"].to_numpy(dtype=np.float64, copy=False))
    area_ok = np.isfinite(pyr_log_area)

    reactive_xy = reactive[["x", "y"]].to_numpy(dtype=np.float32, copy=False)
    tree = cKDTree(reactive_xy)

    for spec in VARIATIONS.values():
        radius = float(spec["radius_px"])
        dists, _ = tree.query(pyr_xy, k=1, distance_upper_bound=radius)
        selected = np.isfinite(dists) & area_ok
        selected_values = pyr_log_area[selected]
        result[spec["selected_count_column"]] = float(selected.sum())
        if selected.sum() < MIN_SELECTED_CELLS:
            result[spec["feature_column"]] = float("nan")
        else:
            result[spec["feature_column"]] = float(np.median(selected_values))

    return result


def compute_donor_score(*, donor_id: str, data_root: str | Path):
    data_root = Path(data_root)
    cohort = _prepare_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    feature_values = _extract_slide_features(slide_path)
    variation_name = _load_best_variation_from_sidecar()
    feature_column = VARIATIONS[variation_name]["feature_column"]
    value = feature_values.get(feature_column, float("nan"))
    if value is None or not np.isfinite(value):
        return None
    return float(value)


def _extract_all_features(cohort: pd.DataFrame, data_root: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for row in cohort.itertuples(index=False):
        slide_path = data_root / str(row.slide_name)
        feature_values = _extract_slide_features(slide_path)
        rows.append(
            {
                "donor_id": row.donor_id,
                "slide_name": row.slide_name,
                **feature_values,
            }
        )
    return pd.DataFrame(rows)


def _loo_single_feature(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
) -> dict[str, object]:
    frame = df.loc[:, ["donor_id", feature_col, outcome_col, *confounds]].dropna().reset_index(drop=True)
    per_donor: list[dict[str, object]] = []
    pred_resid: list[float] = []
    true_resid: list[float] = []

    if len(frame) < 3:
        return {"loo_predictive_r": float("nan"), "per_donor": per_donor}

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train_conf = _design_matrix(train[confounds])
        x_test_conf = _design_matrix(test[confounds])

        yf_train = train[feature_col].to_numpy(dtype=float)
        yo_train = train[outcome_col].to_numpy(dtype=float)

        beta_feature, *_ = np.linalg.lstsq(x_train_conf, yf_train, rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train_conf, yo_train, rcond=None)

        resid_feature_train = yf_train - x_train_conf @ beta_feature
        resid_outcome_train = yo_train - x_train_conf @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            pred_resid_i = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            resid_feature_test = test[feature_col].to_numpy(dtype=float) - x_test_conf @ beta_feature
            pred_resid_i = float(slope * resid_feature_test[0])

        true_resid_i = float(
            test[outcome_col].to_numpy(dtype=float)[0] - (x_test_conf @ beta_outcome)[0]
        )
        raw_pred_i = float((x_test_conf @ beta_outcome)[0] + pred_resid_i) if np.isfinite(pred_resid_i) else float("nan")

        pred_resid.append(pred_resid_i)
        true_resid.append(true_resid_i)
        per_donor.append(
            {
                "donor_id": str(test["donor_id"].iloc[0]),
                "outcome": float(test[outcome_col].iloc[0]),
                "predicted": raw_pred_i,
                feature_col: float(test[feature_col].iloc[0]),
            }
        )

    pred_arr = np.asarray(pred_resid, dtype=float)
    true_arr = np.asarray(true_resid, dtype=float)
    mask = np.isfinite(pred_arr) & np.isfinite(true_arr)
    if mask.sum() < 3 or np.std(pred_arr[mask]) == 0 or np.std(true_arr[mask]) == 0:
        loo_r = float("nan")
    else:
        loo_r = float(np.corrcoef(pred_arr[mask], true_arr[mask])[0, 1])

    return {"loo_predictive_r": loo_r, "per_donor": per_donor}


def _evaluate_variation(
    table: pd.DataFrame,
    *,
    variation_name: str,
    feature_col: str,
) -> dict[str, object]:
    analysis_cols = ["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUNDS]
    clean = table.loc[:, analysis_cols].dropna().reset_index(drop=True)
    stats = partial_correlation(
        clean,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
    )
    partial_r = float(stats["partial_r"]) if "partial_r" in stats else float("nan")
    p_value = float(stats["p_value"]) if "p_value" in stats else float("nan")
    n_total = int(len(table))
    n_analyzable = int(len(clean))
    coverage_ratio = (n_analyzable / n_total) if n_total else float("nan")
    selection_score = abs(partial_r) * coverage_ratio if np.isfinite(partial_r) else float("nan")

    loo = _loo_single_feature(
        clean,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
    )
    loo_r = float(loo["loo_predictive_r"])
    gap = max(0.0, abs(partial_r) - abs(loo_r)) if np.isfinite(partial_r) and np.isfinite(loo_r) else float("nan")
    penalty = max(0.0, gap - 0.20) if np.isfinite(gap) else float("nan")
    adjusted_score = selection_score - penalty if np.isfinite(selection_score) and np.isfinite(penalty) else selection_score

    return {
        "name": variation_name,
        "description": VARIATIONS[variation_name]["description"],
        "radius_px": float(VARIATIONS[variation_name]["radius_px"]),
        "feature_column": feature_col,
        "selected_count_column": VARIATIONS[variation_name]["selected_count_column"],
        "partial_r": partial_r,
        "p_value": p_value,
        "selection_score": selection_score,
        "loo_predictive_r": loo_r,
        "is_loo_gap": gap,
        "penalty": penalty,
        "adjusted_score": adjusted_score,
        "coverage_ratio": coverage_ratio,
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "donor_ids_used": clean["donor_id"].tolist(),
        "per_donor_loo": loo["per_donor"],
    }


def _rank_variations(results: list[dict[str, object]]) -> list[dict[str, object]]:
    return sorted(
        results,
        key=lambda row: (
            -(-1e9 if row["selection_score"] is None or not math.isfinite(float(row["selection_score"])) else float(row["selection_score"])),
            -(-1e9 if row["adjusted_score"] is None or not math.isfinite(float(row["adjusted_score"])) else float(row["adjusted_score"])),
            -(-1e9 if row["loo_predictive_r"] is None or not math.isfinite(float(row["loo_predictive_r"])) else float(row["loo_predictive_r"])),
        ),
    )


def main() -> int:
    data_root = DATA_ROOT_DEFAULT
    cohort = _prepare_cohort(data_root)
    feature_table = _extract_all_features(cohort, data_root)
    donor_table = cohort.merge(feature_table, on=["donor_id", "slide_name"], how="left")

    variation_results: list[dict[str, object]] = []
    for variation_name, spec in VARIATIONS.items():
        variation_results.append(
            _evaluate_variation(
                donor_table,
                variation_name=variation_name,
                feature_col=spec["feature_column"],
            )
        )
    ranked = _rank_variations(variation_results)
    best = ranked[0]

    best_feature_col = str(best["feature_column"])
    donor_feature_path = Path("/scratch/donor_feature_table.csv")
    extra_columns = []
    for spec in VARIATIONS.values():
        if spec["feature_column"] != best_feature_col:
            extra_columns.append(spec["feature_column"])
    extra_columns.extend(KEY_COUNT_COLUMNS)
    write_donor_feature_table(
        donor_feature_path,
        donor_table,
        feature_column=best_feature_col,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUNDS,
        extra_columns=extra_columns,
    )

    results = build_results_payload(
        status="ok",
        feature_name=best_feature_col,
        outcome=OUTCOME_COLUMN,
        n_total=int(best["n_total"]),
        n_analyzable=int(best["n_analyzable"]),
        partial_r=_safe_float(best["partial_r"]),
        ci_lo=None,
        ci_hi=None,
        p_value=_safe_float(best["p_value"]),
        loo_predictive_r=_safe_float(best["loo_predictive_r"]),
        donor_ids_used=list(best["donor_ids_used"]),
        covariates=CONFOUNDS,
        recomputed_from_raw=True,
        registry_written=False,
        artifacts={"donor_feature_table": str(donor_feature_path)},
        best_variation=str(best["name"]),
        ranked_variations=[
            {
                k: (_safe_float(v) if k in {
                    "radius_px",
                    "partial_r",
                    "p_value",
                    "selection_score",
                    "loo_predictive_r",
                    "is_loo_gap",
                    "penalty",
                    "adjusted_score",
                    "coverage_ratio",
                } else v)
                for k, v in row.items()
                if k != "per_donor_loo"
            }
            for row in ranked
        ],
        feature_column=best_feature_col,
        selection_score=_safe_float(best["selection_score"]),
        adjusted_score=_safe_float(best["adjusted_score"]),
        coverage_ratio=_safe_float(best["coverage_ratio"]),
        is_loo_gap=_safe_float(best["is_loo_gap"]),
        gap_penalty=_safe_float(best["penalty"]),
    )
    write_results_payload(Path("/scratch/results.json"), results)

    print(f"HYPOTHESIS FAMILY: {FEATURE_FAMILY}")
    print(f"BEST VARIATION: {best['name']}")
    print(f"  IS partial r:      {_fmt(best['partial_r'])}")
    print(f"  Selection score:   {_fmt(best['selection_score'])}")
    print(f"  LOO predictive r:  {_fmt(best['loo_predictive_r'])}  (diagnostic)")
    print(f"  IS-LOO Gap:        {_fmt(best['is_loo_gap'])}  (penalty={_fmt(best['penalty'])})")
    print(f"  Adjusted Score:    {_fmt(best['adjusted_score'])}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for row in ranked:
        print(
            f"  {row['name']}  {_fmt(row['partial_r'])}  "
            f"{_fmt(row['selection_score'])}  {_fmt(row['loo_predictive_r'])}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best_feature_col}")
    for row in best["per_donor_loo"]:
        print(
            f"  {row['donor_id']}  {_fmt(row['outcome'])}  {_fmt(row['predicted'])}  "
            f"{_fmt(row[best_feature_col])}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
PY

python /scratch/result.py | tee /scratch/run_stdout.txt
