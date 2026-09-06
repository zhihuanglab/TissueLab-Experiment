from __future__ import annotations

import json
import math
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, "/shared/lib")
from shared_analysis.artifacts import write_donor_feature_table
from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort
from shared_analysis.stats import partial_correlation

# Canonical replay target for the evaluator.
HYPOTHESIS_FAMILY = "ca1_reactive_local_purity"
BEST_VARIATION = "reactive_purity_r60um"
FEATURE_NAME = "ca1_reactive_local_purity"
FEATURE_COLUMN = "ca1_reactive_local_purity_r60um"

OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]

ASTRO_LINEAGE = {"Astrocyte", "Reactive Astrocyte"}
REACTIVE_LABEL = "Reactive Astrocyte"
TARGET_REGION = "CA1"

VARIATIONS = {
    "reactive_purity_r40um": {
        "radius_px": 80.0,
        "feature_column": "ca1_reactive_local_purity_r40um",
        "description": "Mean local CA1 reactive purity around each CA1 Reactive Astrocyte within 40 µm (~80 px).",
    },
    "reactive_purity_r60um": {
        "radius_px": 120.0,
        "feature_column": "ca1_reactive_local_purity_r60um",
        "description": "Mean local CA1 reactive purity around each CA1 Reactive Astrocyte within 60 µm (~120 px).",
    },
}


def _suppress_zarr_sidecar_warnings() -> None:
    warnings.filterwarnings(
        "ignore",
        message=r"Object at .* is not recognized as a component of a Zarr hierarchy\.",
        category=UserWarning,
    )


def _design_matrix(df: pd.DataFrame) -> np.ndarray:
    x = np.column_stack([np.ones(len(df), dtype=float)] + [df[col].to_numpy(dtype=float) for col in df.columns])
    return x


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


def _json_default(obj):
    if isinstance(obj, (np.floating,)):
        val = float(obj)
        return None if not math.isfinite(val) else val
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Not JSON serializable: {type(obj)!r}")


def _local_reactive_purity_from_cells(cells: pd.DataFrame, radius_px: float) -> tuple[float, dict[str, int | float | None]]:
    ca1 = cells.loc[(cells["region"] == TARGET_REGION) & (cells["cell_type"].isin(ASTRO_LINEAGE)), ["x", "y", "cell_type"]].copy()
    if ca1.empty:
        return float("nan"), {"ca1_astro_lineage_n": 0, "ca1_reactive_n": 0, "valid_centers": 0}

    coords = ca1[["x", "y"]].to_numpy(dtype=np.float32)
    reactive_mask = (ca1["cell_type"].to_numpy() == REACTIVE_LABEL)
    reactive_indices = np.flatnonzero(reactive_mask)
    reactive_n = int(reactive_indices.size)
    if reactive_n == 0:
        return float("nan"), {
            "ca1_astro_lineage_n": int(len(ca1)),
            "ca1_reactive_n": 0,
            "valid_centers": 0,
        }

    tree = cKDTree(coords)
    neighbor_lists = tree.query_ball_point(coords[reactive_indices], r=radius_px)

    purities: list[float] = []
    for center_idx, neighbors in zip(reactive_indices, neighbor_lists):
        total_neighbors = len(neighbors)
        if center_idx in neighbors:
            total_neighbors -= 1
        if total_neighbors < 2:
            continue
        reactive_neighbors = int(reactive_mask[np.asarray(neighbors, dtype=int)].sum())
        if center_idx in neighbors:
            reactive_neighbors -= 1
        purities.append(reactive_neighbors / total_neighbors)

    if not purities:
        return float("nan"), {
            "ca1_astro_lineage_n": int(len(ca1)),
            "ca1_reactive_n": reactive_n,
            "valid_centers": 0,
        }

    return float(np.mean(purities)), {
        "ca1_astro_lineage_n": int(len(ca1)),
        "ca1_reactive_n": reactive_n,
        "valid_centers": int(len(purities)),
    }


def extract_donor_features(*, slide_path: str | Path) -> dict[str, float | int | None]:
    _suppress_zarr_sidecar_warnings()
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    result: dict[str, float | int | None] = {}
    shared_counts_written = False
    for variation_name, spec in VARIATIONS.items():
        value, counts = _local_reactive_purity_from_cells(cells, radius_px=float(spec["radius_px"]))
        result[spec["feature_column"]] = value
        if not shared_counts_written:
            result["ca1_astro_lineage_n"] = counts["ca1_astro_lineage_n"]
            result["ca1_reactive_n"] = counts["ca1_reactive_n"]
            shared_counts_written = True
        result[f"{variation_name}_valid_centers"] = counts["valid_centers"]
    return result


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Return the canonical winning biomarker score for one donor.

    This replay target recomputes from the raw .zarr slide under data_root.
    """
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    features = extract_donor_features(slide_path=slide_path)
    value = features.get(FEATURE_COLUMN)
    return None if value is None or not math.isfinite(float(value)) else float(value)


@dataclass
class VariationEvaluation:
    name: str
    feature_column: str
    partial_r: float | None
    p_value: float | None
    n_analyzable: int
    n_total: int
    selection_score: float | None
    loo_predictive_r: float | None
    is_loo_gap: float | None
    penalty: float | None
    adjusted_score: float | None
    loo_rows: list[dict[str, object]]


def prepare_analysis_table(data_root: Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map({"Female": 0.0, "Male": 1.0}).astype(float)
    rows = []
    for row in cohort.itertuples(index=False):
        slide_name = str(row.slide_name)
        donor_id = str(row.donor_id)
        slide_path = data_root / slide_name
        features = extract_donor_features(slide_path=slide_path)
        rows.append({"donor_id": donor_id, "slide_name": slide_name, **features})
    feat_df = pd.DataFrame(rows)
    merged = cohort.merge(feat_df, on=["donor_id", "slide_name"], how="left")
    return merged


def loo_predictions(df: pd.DataFrame, *, feature_col: str) -> tuple[float | None, list[dict[str, object]]]:
    needed = ["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]
    frame = df.loc[:, needed].dropna().reset_index(drop=True)
    if len(frame) < 3:
        return None, []

    preds_resid: list[float] = []
    actuals_resid: list[float] = []
    rows: list[dict[str, object]] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train_conf = _design_matrix(train[CONFOUND_COLUMNS])
        x_test_conf = _design_matrix(test[CONFOUND_COLUMNS])

        beta_feature, *_ = np.linalg.lstsq(x_train_conf, train[feature_col].to_numpy(dtype=float), rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train_conf, train[OUTCOME_COLUMN].to_numpy(dtype=float), rcond=None)

        resid_feature_train = train[feature_col].to_numpy(dtype=float) - x_train_conf @ beta_feature
        resid_outcome_train = train[OUTCOME_COLUMN].to_numpy(dtype=float) - x_train_conf @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            pred_resid = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            resid_feature_test = test[feature_col].to_numpy(dtype=float) - x_test_conf @ beta_feature
            pred_resid = float(slope * resid_feature_test[0])

        actual_resid = float(test[OUTCOME_COLUMN].to_numpy(dtype=float)[0] - (x_test_conf @ beta_outcome)[0])
        pred_outcome = float((x_test_conf @ beta_outcome)[0] + pred_resid) if math.isfinite(pred_resid) else float("nan")

        preds_resid.append(pred_resid)
        actuals_resid.append(actual_resid)
        rows.append(
            {
                "donor_id": str(test.loc[0, "donor_id"]),
                "outcome": float(test.loc[0, OUTCOME_COLUMN]),
                "predicted": pred_outcome,
                "predicted_residual": pred_resid,
                "actual_residual": actual_resid,
                feature_col: float(test.loc[0, feature_col]),
            }
        )

    pred_arr = np.asarray(preds_resid, dtype=float)
    act_arr = np.asarray(actuals_resid, dtype=float)
    mask = np.isfinite(pred_arr) & np.isfinite(act_arr)
    if mask.sum() < 3 or np.std(pred_arr[mask]) == 0 or np.std(act_arr[mask]) == 0:
        loo_r = None
    else:
        loo_r = float(np.corrcoef(pred_arr[mask], act_arr[mask])[0, 1])
    return loo_r, rows


def evaluate_variation(df: pd.DataFrame, *, variation_name: str, feature_col: str) -> VariationEvaluation:
    clean = df.loc[:, ["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]].dropna().reset_index(drop=True)
    n_total = int(len(df))
    n_analyzable = int(len(clean))

    if n_analyzable < 3:
        return VariationEvaluation(
            name=variation_name,
            feature_column=feature_col,
            partial_r=None,
            p_value=None,
            n_analyzable=n_analyzable,
            n_total=n_total,
            selection_score=None,
            loo_predictive_r=None,
            is_loo_gap=None,
            penalty=None,
            adjusted_score=None,
            loo_rows=[],
        )

    pc = partial_correlation(
        clean,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    partial_r = _safe_float(pc.get("partial_r"))
    p_value = _safe_float(pc.get("p_value"))
    selection_score = None if partial_r is None else abs(partial_r) * (n_analyzable / n_total)

    loo_r, loo_rows = loo_predictions(clean, feature_col=feature_col)
    gap = None
    penalty = None
    adjusted = None
    if partial_r is not None and loo_r is not None:
        gap = abs(partial_r) - abs(loo_r)
        penalty = max(0.0, gap - 0.15) / 2.0
        adjusted = abs(loo_r) - penalty

    return VariationEvaluation(
        name=variation_name,
        feature_column=feature_col,
        partial_r=partial_r,
        p_value=p_value,
        n_analyzable=n_analyzable,
        n_total=n_total,
        selection_score=selection_score,
        loo_predictive_r=loo_r,
        is_loo_gap=gap,
        penalty=penalty,
        adjusted_score=adjusted,
        loo_rows=loo_rows,
    )


def ranked_variations_payload(evals: Iterable[VariationEvaluation]) -> list[dict[str, object]]:
    rows = []
    for ev in evals:
        rows.append(
            {
                "name": ev.name,
                "feature_column": ev.feature_column,
                "partial_r": ev.partial_r,
                "selection_score": ev.selection_score,
                "loo_predictive_r": ev.loo_predictive_r,
                "is_loo_gap": ev.is_loo_gap,
                "gap": ev.is_loo_gap,
                "penalty": ev.penalty,
                "gap_penalty": ev.penalty,
                "adjusted_score": ev.adjusted_score,
                "p_value": ev.p_value,
                "n_analyzable": ev.n_analyzable,
                "n_total": ev.n_total,
            }
        )
    rows.sort(
        key=lambda d: (
            -1e9 if d["selection_score"] is None else -float(d["selection_score"]),
            -1e9 if d["adjusted_score"] is None else -float(d["adjusted_score"]),
            -1e9 if d["partial_r"] is None else -abs(float(d["partial_r"])),
            d["name"],
        )
    )
    return rows


def fmt(x: float | None) -> str:
    if x is None or not math.isfinite(float(x)):
        return "nan"
    return f"{float(x):.4f}"


def main() -> None:
    data_root = Path("/data")
    analysis_table = prepare_analysis_table(data_root)
    donor_feature_table_path = Path("/scratch/donor_feature_table.csv")
    extra_cols = [
        "ca1_astro_lineage_n",
        "ca1_reactive_n",
        *[spec["feature_column"] for spec in VARIATIONS.values()],
        *[f"{name}_valid_centers" for name in VARIATIONS],
    ]
    write_donor_feature_table(
        donor_feature_table_path,
        analysis_table,
        feature_column=FEATURE_COLUMN,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUND_COLUMNS,
        extra_columns=[c for c in extra_cols if c != FEATURE_COLUMN],
    )

    evaluations: list[VariationEvaluation] = []
    for variation_name, spec in VARIATIONS.items():
        evaluations.append(evaluate_variation(analysis_table, variation_name=variation_name, feature_col=spec["feature_column"]))

    ranked = ranked_variations_payload(evaluations)
    best_name = ranked[0]["name"]
    best_eval = next(ev for ev in evaluations if ev.name == best_name)

    payload = {
        "status": "ok",
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "feature_name": FEATURE_NAME,
        "best_variation": best_name,
        "feature_column": best_eval.feature_column,
        "outcome": OUTCOME_COLUMN,
        "n_total": best_eval.n_total,
        "n_analyzable": best_eval.n_analyzable,
        "partial_r": best_eval.partial_r,
        "selection_score": best_eval.selection_score,
        "loo_predictive_r": best_eval.loo_predictive_r,
        "is_loo_gap": best_eval.is_loo_gap,
        "penalty": best_eval.penalty,
        "adjusted_score": best_eval.adjusted_score,
        "p_value": best_eval.p_value,
        "covariates": CONFOUND_COLUMNS,
        "donor_ids_used": [
            str(x)
            for x in analysis_table.loc[analysis_table[best_eval.feature_column].notna(), "donor_id"].tolist()
        ],
        "ranked_variations": ranked,
        "donor_feature_table": str(donor_feature_table_path),
        "loo_rows": best_eval.loo_rows,
    }
    with open("/scratch/results.json", "w") as f:
        json.dump(payload, f, indent=2, default=_json_default)

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best_name}")
    print(f"  IS partial r:      {fmt(best_eval.partial_r)}")
    print(f"  Selection score:   {fmt(best_eval.selection_score)}")
    print(f"  LOO predictive r:  {fmt(best_eval.loo_predictive_r)}  (diagnostic)")
    print(f"  IS-LOO Gap:        {fmt(best_eval.is_loo_gap)}  (penalty={fmt(best_eval.penalty)})")
    print(f"  Adjusted Score:    {fmt(best_eval.adjusted_score)}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for row in ranked:
        print(
            f"  {row['name']}  {fmt(row['partial_r'])}  {fmt(row['selection_score'])}  {fmt(row['loo_predictive_r'])}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best_eval.feature_column}")
    for row in best_eval.loo_rows:
        pred = row['predicted']
        pred_text = f"{float(pred):.4f}" if pred is not None and math.isfinite(float(pred)) else "nan"
        print(
            f"  {row['donor_id']}  {row['outcome']:.4f}  "
            f"{pred_text}  "
            f"{row[best_eval.feature_column]:.6f}"
        )


if __name__ == "__main__":
    main()
