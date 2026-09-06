from __future__ import annotations

import json
import math
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, "/shared/lib")
from shared_analysis import build_cell_table, load_training_cohort  # noqa: E402
from shared_analysis.artifacts import write_donor_feature_table  # noqa: E402
from shared_analysis.stats import partial_correlation  # noqa: E402


DATA_ROOT = Path("/data")
OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
HYPOTHESIS_FAMILY = "ca1_pyramidal_reactive_astrocyte_neighbor_fraction"
VARIATION_RADII = {
    "r4_neighbor_fraction": 4.0,
    "r6_neighbor_fraction": 6.0,
}
FEATURE_COLUMNS = {
    name: f"ca1_pyramidal_reactive_neighbor_fraction_{name}"
    for name in VARIATION_RADII
}
CANONICAL_VARIATION = "r6_neighbor_fraction"


warnings.filterwarnings(
    "ignore",
    message=r"Object at .* is not recognized as a component of a Zarr hierarchy.",
)
warnings.filterwarnings("ignore", category=RuntimeWarning)


def _safe_float(value) -> float:
    try:
        val = float(value)
    except Exception:
        return float("nan")
    if math.isnan(val) or math.isinf(val):
        return float("nan")
    return val


def _load_cohort(data_root: str | Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map({"Female": 0.0, "Male": 1.0})
    return cohort


def _extract_slide_metrics(slide_path: str | Path) -> dict[str, float]:
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    ca1 = cells.loc[cells["region"] == "CA1", ["x", "y", "cell_type"]].copy()
    pyramidal = ca1.loc[ca1["cell_type"] == "Pyramidal Neuron", ["x", "y"]].to_numpy(dtype=float)
    reactive = ca1.loc[ca1["cell_type"] == "Reactive Astrocyte", ["x", "y"]].to_numpy(dtype=float)

    result: dict[str, float] = {
        "ca1_cell_count": float(len(ca1)),
        "ca1_pyramidal_count": float(len(pyramidal)),
        "ca1_reactive_astrocyte_count": float(len(reactive)),
    }

    if len(pyramidal) == 0:
        nearest_d = np.asarray([], dtype=float)
        for variation in VARIATION_RADII:
            result[FEATURE_COLUMNS[variation]] = float("nan")
    elif len(reactive) == 0:
        nearest_d = np.full(len(pyramidal), np.inf, dtype=float)
        for variation in VARIATION_RADII:
            result[FEATURE_COLUMNS[variation]] = 0.0
    else:
        nearest_d, _ = cKDTree(reactive).query(pyramidal, k=1)
        nearest_d = np.asarray(nearest_d, dtype=float)
        for variation, radius in VARIATION_RADII.items():
            result[FEATURE_COLUMNS[variation]] = float(np.mean(nearest_d <= radius))

    result["ca1_pyramidal_nearest_reactive_min_distance"] = (
        float(np.nanmin(nearest_d)) if nearest_d.size else float("nan")
    )
    result["ca1_pyramidal_nearest_reactive_median_distance"] = (
        float(np.nanmedian(nearest_d)) if nearest_d.size else float("nan")
    )
    return result


def _resolve_canonical_variation() -> str:
    if CANONICAL_VARIATION in VARIATION_RADII:
        return CANONICAL_VARIATION
    results_path = Path(__file__).with_name("results.json")
    if results_path.exists():
        try:
            payload = json.loads(results_path.read_text())
            variation = payload.get("best_variation")
            if variation in VARIATION_RADII:
                return str(variation)
        except Exception:
            pass
    return "r4_neighbor_fraction"


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Return the winning CA1 pyramidal reactive-neighbor fraction for one donor.

    The canonical variation is rewritten after the sweep runs so held-out
    evaluation can replay the local winner from raw data.
    """
    data_root = Path(data_root)
    cohort = _load_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    slide_metrics = _extract_slide_metrics(slide_path)
    feature_col = FEATURE_COLUMNS[_resolve_canonical_variation()]
    value = slide_metrics.get(feature_col, float("nan"))
    return None if not np.isfinite(value) else float(value)


def _design_matrix(frame: pd.DataFrame, *, confounds: list[str], feature_col: str | None = None) -> np.ndarray:
    cols = [np.ones(len(frame), dtype=float)]
    for col in confounds:
        cols.append(frame[col].to_numpy(dtype=float))
    if feature_col is not None:
        cols.append(frame[feature_col].to_numpy(dtype=float))
    return np.column_stack(cols)


def _safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return float("nan")
    a = a[mask]
    b = b[mask]
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _loo_for_feature(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
) -> tuple[float, list[dict[str, float]]]:
    required = ["donor_id", outcome_col, feature_col, *confounds]
    frame = df.loc[:, [c for c in df.columns if c in set(required) or c.startswith("ca1_")]].copy()
    frame = frame.dropna(subset=required).reset_index(drop=True)
    records: list[dict[str, float]] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train_conf = _design_matrix(train, confounds=confounds)
        x_test_conf = _design_matrix(test, confounds=confounds)

        beta_feature, *_ = np.linalg.lstsq(
            x_train_conf,
            train[feature_col].to_numpy(dtype=float),
            rcond=None,
        )
        beta_outcome_conf, *_ = np.linalg.lstsq(
            x_train_conf,
            train[outcome_col].to_numpy(dtype=float),
            rcond=None,
        )

        resid_feature_train = train[feature_col].to_numpy(dtype=float) - x_train_conf @ beta_feature
        resid_outcome_train = train[outcome_col].to_numpy(dtype=float) - x_train_conf @ beta_outcome_conf
        denom = float(np.dot(resid_feature_train, resid_feature_train))
        slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom) if denom > 0 else 0.0

        resid_feature_test = test[feature_col].to_numpy(dtype=float) - x_test_conf @ beta_feature
        resid_outcome_test = test[outcome_col].to_numpy(dtype=float) - x_test_conf @ beta_outcome_conf

        residual_pred = float(slope * resid_feature_test[0])
        raw_pred = float((x_test_conf @ beta_outcome_conf)[0] + residual_pred)

        row = test.iloc[0]
        records.append(
            {
                "donor_id": str(row["donor_id"]),
                "outcome": float(row[outcome_col]),
                "predicted": raw_pred,
                "residualized_outcome": float(resid_outcome_test[0]),
                "residualized_predicted": residual_pred,
                "feature_value": float(row[feature_col]),
                "ca1_cell_count": _safe_float(row.get("ca1_cell_count", float("nan"))),
                "ca1_pyramidal_count": _safe_float(row.get("ca1_pyramidal_count", float("nan"))),
                "ca1_reactive_astrocyte_count": _safe_float(row.get("ca1_reactive_astrocyte_count", float("nan"))),
                "ca1_pyramidal_nearest_reactive_min_distance": _safe_float(
                    row.get("ca1_pyramidal_nearest_reactive_min_distance", float("nan"))
                ),
                "ca1_pyramidal_nearest_reactive_median_distance": _safe_float(
                    row.get("ca1_pyramidal_nearest_reactive_median_distance", float("nan"))
                ),
            }
        )

    loo_r = _safe_corr(
        np.asarray([r["residualized_predicted"] for r in records], dtype=float),
        np.asarray([r["residualized_outcome"] for r in records], dtype=float),
    )
    return loo_r, records


def _gap_penalty(is_partial_r: float, loo_r: float) -> tuple[float, float]:
    if not np.isfinite(is_partial_r) or not np.isfinite(loo_r):
        return float("nan"), 0.0
    gap = float(abs(abs(is_partial_r) - abs(loo_r)))
    penalty = max(0.0, gap - 0.15) * 0.5
    return gap, penalty


def _evaluate_variation(
    donor_table: pd.DataFrame,
    *,
    variation: str,
    outcome_col: str,
    confounds: list[str],
) -> dict[str, object]:
    feature_col = FEATURE_COLUMNS[variation]
    required = [feature_col, outcome_col, *confounds]
    frame = donor_table.dropna(subset=required).copy()

    is_stats = partial_correlation(
        donor_table,
        feature_col=feature_col,
        outcome_col=outcome_col,
        confounds=confounds,
    )
    partial_r = _safe_float(is_stats.get("partial_r", float("nan")))
    p_value = _safe_float(is_stats.get("p_value", float("nan")))
    n_analyzable = int(frame.shape[0])
    n_total = int(donor_table.shape[0])
    coverage = n_analyzable / n_total if n_total else float("nan")
    selection_score = abs(partial_r) * coverage if np.isfinite(partial_r) else float("nan")
    loo_r, loo_records = _loo_for_feature(
        donor_table,
        feature_col=feature_col,
        outcome_col=outcome_col,
        confounds=confounds,
    )
    gap, penalty = _gap_penalty(partial_r, loo_r)
    adjusted_score = selection_score - penalty if np.isfinite(selection_score) else float("nan")

    return {
        "variation_name": variation,
        "feature_column": feature_col,
        "radius": VARIATION_RADII[variation],
        "partial_r": partial_r,
        "p_value": p_value,
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "coverage": coverage,
        "selection_score": selection_score,
        "loo_predictive_r": loo_r,
        "is_loo_gap": gap,
        "gap_penalty": penalty,
        "adjusted_score": adjusted_score,
        "loo_records": loo_records,
    }


def _update_script_canonical_variation(script_path: Path, variation: str) -> None:
    text = script_path.read_text()
    new_text, count = re.subn(
        r'CANONICAL_VARIATION = ".*?"',
        f'CANONICAL_VARIATION = "{variation}"',
        text,
        count=1,
    )
    if count == 1 and new_text != text:
        script_path.write_text(new_text)


def _json_ready(obj):
    if isinstance(obj, dict):
        return {str(k): _json_ready(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_json_ready(v) for v in obj]
    if isinstance(obj, (np.floating,)):
        obj = float(obj)
    if isinstance(obj, (np.integer,)):
        obj = int(obj)
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
        return obj
    return obj


def _format_metric(value: float) -> str:
    return "nan" if not np.isfinite(value) else f"{value:.4f}"


def main() -> None:
    cohort = _load_cohort(DATA_ROOT)
    donor_rows: list[dict[str, object]] = []
    for row in cohort.itertuples(index=False):
        slide_metrics = _extract_slide_metrics(DATA_ROOT / str(row.slide_name))
        base = {
            "donor_id": str(row.donor_id),
            "slide_name": str(row.slide_name),
            OUTCOME_COLUMN: _safe_float(getattr(row, OUTCOME_COLUMN)),
            "max_age_vis": _safe_float(row.max_age_vis),
            "braak_numeric": _safe_float(row.braak_numeric),
            "cerad_ordinal": _safe_float(row.cerad_ordinal),
            "sex_binary": _safe_float({"Female": 0.0, "Male": 1.0}.get(row.sex, np.nan)),
        }
        base.update(slide_metrics)
        donor_rows.append(base)

    donor_table = pd.DataFrame(donor_rows)

    evaluations = [
        _evaluate_variation(
            donor_table,
            variation=variation,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUND_COLUMNS,
        )
        for variation in VARIATION_RADII
    ]
    ranked = sorted(
        evaluations,
        key=lambda d: (
            -np.inf if not np.isfinite(d["selection_score"]) else -float(d["selection_score"]),
            -np.inf if not np.isfinite(d["adjusted_score"]) else -float(d["adjusted_score"]),
            str(d["variation_name"]),
        ),
    )
    if not ranked:
        raise RuntimeError("No variations were evaluated.")
    best = ranked[0]
    best_variation = str(best["variation_name"])
    feature_column = str(best["feature_column"])

    _update_script_canonical_variation(Path(__file__), best_variation)

    donor_feature_table_path = Path("/scratch/donor_feature_table.csv")
    extra_cols = [
        "ca1_cell_count",
        "ca1_pyramidal_count",
        "ca1_reactive_astrocyte_count",
        "ca1_pyramidal_nearest_reactive_min_distance",
        "ca1_pyramidal_nearest_reactive_median_distance",
        *[FEATURE_COLUMNS[name] for name in VARIATION_RADII if FEATURE_COLUMNS[name] != feature_column],
    ]
    write_donor_feature_table(
        donor_feature_table_path,
        donor_table,
        feature_column=feature_column,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUND_COLUMNS,
        extra_columns=extra_cols,
    )

    payload = {
        "status": "ok",
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": best_variation,
        "feature_name": HYPOTHESIS_FAMILY,
        "feature_column": feature_column,
        "outcome": OUTCOME_COLUMN,
        "covariates": CONFOUND_COLUMNS,
        "n_total": int(best["n_total"]),
        "n_analyzable": int(best["n_analyzable"]),
        "partial_r": _json_ready(best["partial_r"]),
        "p_value": _json_ready(best["p_value"]),
        "selection_score": _json_ready(best["selection_score"]),
        "loo_predictive_r": _json_ready(best["loo_predictive_r"]),
        "is_loo_gap": _json_ready(best["is_loo_gap"]),
        "gap_penalty": _json_ready(best["gap_penalty"]),
        "adjusted_score": _json_ready(best["adjusted_score"]),
        "donor_ids_used": donor_table.dropna(
            subset=[feature_column, OUTCOME_COLUMN, *CONFOUND_COLUMNS]
        )["donor_id"].astype(str).tolist(),
        "artifacts": {
            "donor_feature_table": str(donor_feature_table_path),
        },
        "ranked_variations": [
            {
                k: _json_ready(v)
                for k, v in variation.items()
                if k != "loo_records"
            }
            for variation in ranked
        ],
        "per_donor_loo": _json_ready(best["loo_records"]),
    }
    Path("/scratch/results.json").write_text(json.dumps(_json_ready(payload), indent=2))

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best_variation}")
    print(f"  IS partial r:      {_format_metric(float(best['partial_r']))}")
    print(f"  Selection score:   {_format_metric(float(best['selection_score']))}")
    print(f"  LOO predictive r:  {_format_metric(float(best['loo_predictive_r']))}  (diagnostic)")
    print(
        f"  IS-LOO Gap:        {_format_metric(float(best['is_loo_gap']))}  "
        f"(penalty={_format_metric(float(best['gap_penalty']))})"
    )
    print(f"  Adjusted Score:    {_format_metric(float(best['adjusted_score']))}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for variation in ranked:
        print(
            "  "
            f"{variation['variation_name']}  "
            f"{_format_metric(float(variation['partial_r']))}  "
            f"{_format_metric(float(variation['selection_score']))}  "
            f"{_format_metric(float(variation['loo_predictive_r']))}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {feature_column}  ca1_pyramidal_count  ca1_reactive_astrocyte_count")
    for rec in best["loo_records"]:
        print(
            "  "
            f"{rec['donor_id']}  "
            f"{rec['outcome']:.6f}  "
            f"{rec['predicted']:.6f}  "
            f"{rec['feature_value']:.6f}  "
            f"{int(rec['ca1_pyramidal_count']) if np.isfinite(rec['ca1_pyramidal_count']) else 'nan'}  "
            f"{int(rec['ca1_reactive_astrocyte_count']) if np.isfinite(rec['ca1_reactive_astrocyte_count']) else 'nan'}"
        )


if __name__ == "__main__":
    main()
