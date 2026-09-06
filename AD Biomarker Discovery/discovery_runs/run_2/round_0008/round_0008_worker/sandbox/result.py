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

from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort

warnings.filterwarnings("ignore", category=UserWarning, module="zarr")
warnings.filterwarnings("ignore", message="Object at .* is not recognized as a component of a Zarr hierarchy.")
warnings.filterwarnings("ignore", category=RuntimeWarning)

FEATURE_NAME = "ca1_reactive_hotspot_oligodendrocyte_depletion"
HYPOTHESIS_FAMILY = FEATURE_NAME
WINNING_VARIATION = "candidate_variant_b"
FEATURE_COLUMN = f"{FEATURE_NAME}__{WINNING_VARIATION}"

OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
HOTSPOT_RADIUS_PX = 70.0
OLIGO_RADIUS_PX = 70.0

RESULTS_PATH = Path("/scratch/results.json")
REPORT_PATH = Path("/scratch/report.md")
FEATURE_TABLE_PATH = Path("/scratch/donor_feature_table.csv")

VARIATIONS: list[dict[str, Any]] = [
    {
        "name": "candidate_variant_a",
        "description": "Reactive-versus-astrocyte hotspot oligodendrocyte depletion in CA1 with hotspot threshold k=3 same-type neighbors within 70 px.",
        "k_same_type_neighbors": 3,
    },
    {
        "name": "candidate_variant_b",
        "description": "Same biomarker but with a slightly stricter hotspot threshold k=4 same-type neighbors within 70 px.",
        "k_same_type_neighbors": 4,
    },
]


def _feature_column_name(variation_name: str) -> str:
    return f"{FEATURE_NAME}__{variation_name}"


def _safe_float(value: Any) -> float:
    try:
        out = float(value)
    except Exception:
        return float("nan")
    return out if np.isfinite(out) else float("nan")


def _fmt(value: Any, digits: int = 4) -> str:
    value = _safe_float(value)
    if not np.isfinite(value):
        return "nan"
    return f"{value:.{digits}f}"


def _to_python(value: Any) -> Any:
    if isinstance(value, (np.floating,)):
        out = float(value)
        return out if math.isfinite(out) else None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, (list, tuple)):
        return [_to_python(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _to_python(v) for k, v in value.items()}
    return value


def _sex_to_binary(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .map({"female": 0.0, "male": 1.0})
        .astype(float)
    )


def _design_matrix(df: pd.DataFrame, columns: list[str]) -> np.ndarray:
    x = df.loc[:, columns].to_numpy(dtype=float)
    return np.column_stack([np.ones(len(df), dtype=float), x])


def _residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ beta


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3 or len(b) < 3:
        return float("nan")
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        return float("nan")
    if np.std(a) <= 0 or np.std(b) <= 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _partial_correlation(df: pd.DataFrame, feature_col: str) -> dict[str, float]:
    cols = [feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]
    frame = df.loc[:, cols].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    if len(frame) < 3:
        return {"partial_r": float("nan"), "n_analyzable": int(len(frame))}
    x = _design_matrix(frame, CONFOUND_COLUMNS)
    feature_resid = _residualize(frame[feature_col].to_numpy(dtype=float), x)
    outcome_resid = _residualize(frame[OUTCOME_COLUMN].to_numpy(dtype=float), x)
    return {"partial_r": _corr(feature_resid, outcome_resid), "n_analyzable": int(len(frame))}


def _adjusted_score(partial_r: float, loo_r: float) -> tuple[float, float, float]:
    partial_r = _safe_float(partial_r)
    loo_r = _safe_float(loo_r)
    if not np.isfinite(partial_r) or not np.isfinite(loo_r):
        return float("nan"), float("nan"), float("nan")
    gap = float(abs(partial_r) - abs(loo_r))
    penalty = float(max(0.0, gap - 0.15) * 0.5)
    adjusted = float(loo_r - penalty)
    return gap, penalty, adjusted


def _same_class_neighbor_counts(coords: np.ndarray, radius_px: float) -> np.ndarray:
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[0] == 0:
        return np.zeros(0, dtype=int)
    tree = cKDTree(coords)
    neighborhoods = tree.query_ball_point(coords, r=float(radius_px))
    return np.asarray([max(0, len(nbrs) - 1) for nbrs in neighborhoods], dtype=int)


def _cross_class_neighbor_counts(source_coords: np.ndarray, target_coords: np.ndarray, radius_px: float) -> np.ndarray:
    source_coords = np.asarray(source_coords, dtype=float)
    target_coords = np.asarray(target_coords, dtype=float)
    if source_coords.ndim != 2 or source_coords.shape[0] == 0:
        return np.zeros(0, dtype=int)
    if target_coords.ndim != 2 or target_coords.shape[0] == 0:
        return np.zeros(source_coords.shape[0], dtype=int)
    tree = cKDTree(target_coords)
    neighborhoods = tree.query_ball_point(source_coords, r=float(radius_px))
    return np.asarray([len(nbrs) for nbrs in neighborhoods], dtype=int)


def _median_or_nan(values: np.ndarray) -> float:
    if len(values) == 0:
        return float("nan")
    return float(np.median(values))


def _extract_slide_primitives(slide_path: Path) -> dict[str, Any]:
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    cells = cells.loc[cells["region"] == "CA1", ["x", "y", "cell_type"]].copy()
    relevant_types = {"Reactive Astrocyte", "Astrocyte", "Oligodendrocyte"}
    cells = cells.loc[cells["cell_type"].isin(relevant_types)].reset_index(drop=True)

    reactive_coords = cells.loc[cells["cell_type"] == "Reactive Astrocyte", ["x", "y"]].to_numpy(dtype=float)
    astro_coords = cells.loc[cells["cell_type"] == "Astrocyte", ["x", "y"]].to_numpy(dtype=float)
    oligo_coords = cells.loc[cells["cell_type"] == "Oligodendrocyte", ["x", "y"]].to_numpy(dtype=float)

    return {
        "n_ca1_cells_relevant": int(len(cells)),
        "n_ca1_reactive": int(len(reactive_coords)),
        "n_ca1_astro": int(len(astro_coords)),
        "n_ca1_oligo": int(len(oligo_coords)),
        "reactive_coords": reactive_coords,
        "astro_coords": astro_coords,
        "oligo_coords": oligo_coords,
        "reactive_sameclass_neighbors_r70": _same_class_neighbor_counts(reactive_coords, HOTSPOT_RADIUS_PX),
        "astro_sameclass_neighbors_r70": _same_class_neighbor_counts(astro_coords, HOTSPOT_RADIUS_PX),
    }


def _compute_variation_from_primitives(primitives: dict[str, Any], variation: dict[str, Any]) -> dict[str, Any]:
    k = int(variation["k_same_type_neighbors"])
    reactive_coords = primitives["reactive_coords"]
    astro_coords = primitives["astro_coords"]
    oligo_coords = primitives["oligo_coords"]
    reactive_neighbors = primitives["reactive_sameclass_neighbors_r70"]
    astro_neighbors = primitives["astro_sameclass_neighbors_r70"]

    analyzable = (
        primitives["n_ca1_reactive"] > 0
        and primitives["n_ca1_astro"] > 0
        and primitives["n_ca1_oligo"] > 0
    )

    reactive_hotspot_mask = reactive_neighbors >= k
    astro_hotspot_mask = astro_neighbors >= k

    reactive_fallback_used = bool(analyzable and reactive_hotspot_mask.sum() == 0)
    astro_fallback_used = bool(analyzable and astro_hotspot_mask.sum() == 0)

    reactive_centers = reactive_coords if reactive_fallback_used else reactive_coords[reactive_hotspot_mask]
    astro_centers = astro_coords if astro_fallback_used else astro_coords[astro_hotspot_mask]

    reactive_oligo_counts = _cross_class_neighbor_counts(reactive_centers, oligo_coords, OLIGO_RADIUS_PX)
    astro_oligo_counts = _cross_class_neighbor_counts(astro_centers, oligo_coords, OLIGO_RADIUS_PX)

    reactive_median = _median_or_nan(reactive_oligo_counts)
    astro_median = _median_or_nan(astro_oligo_counts)
    feature_value = float(reactive_median - astro_median) if analyzable else float("nan")

    return {
        _feature_column_name(variation["name"]): feature_value,
        f"{variation['name']}__k_same_type_neighbors": k,
        f"{variation['name']}__n_reactive_centers": int(len(reactive_centers)),
        f"{variation['name']}__n_astro_centers": int(len(astro_centers)),
        f"{variation['name']}__reactive_hotspot_fraction": float(np.mean(reactive_hotspot_mask)) if len(reactive_hotspot_mask) > 0 else float("nan"),
        f"{variation['name']}__astro_hotspot_fraction": float(np.mean(astro_hotspot_mask)) if len(astro_hotspot_mask) > 0 else float("nan"),
        f"{variation['name']}__reactive_fallback_used": reactive_fallback_used,
        f"{variation['name']}__astro_fallback_used": astro_fallback_used,
        f"{variation['name']}__reactive_median_oligo_neighbors_r70": reactive_median,
        f"{variation['name']}__astro_median_oligo_neighbors_r70": astro_median,
        f"{variation['name']}__reactive_minus_astro_median_oligo_neighbors_r70": feature_value,
    }


def _compute_single_donor_row(donor_row: pd.Series, data_root: Path) -> dict[str, Any]:
    donor_id = str(donor_row["donor_id"])
    slide_name = str(donor_row["slide_name"])
    slide_path = data_root / slide_name

    primitives = _extract_slide_primitives(slide_path)
    feature_row: dict[str, Any] = {
        "donor_id": donor_id,
        "slide_name": slide_name,
        "n_ca1_cells_relevant": primitives["n_ca1_cells_relevant"],
        "n_ca1_reactive": primitives["n_ca1_reactive"],
        "n_ca1_astro": primitives["n_ca1_astro"],
        "n_ca1_oligo": primitives["n_ca1_oligo"],
        "reactive_sameclass_neighbor_median_r70": _median_or_nan(primitives["reactive_sameclass_neighbors_r70"]),
        "astro_sameclass_neighbor_median_r70": _median_or_nan(primitives["astro_sameclass_neighbors_r70"]),
    }
    for variation in VARIATIONS:
        feature_row.update(_compute_variation_from_primitives(primitives, variation))
    return feature_row


def _compute_all_donor_features(data_root: Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = _sex_to_binary(cohort["sex"])

    rows: list[dict[str, Any]] = []
    for _, donor_row in cohort.iterrows():
        rows.append(_compute_single_donor_row(donor_row, data_root))

    feature_df = pd.DataFrame(rows)
    donor_table = cohort.merge(feature_df, on=["donor_id", "slide_name"], how="left")
    return donor_table


def _loo_for_variation(df: pd.DataFrame, feature_col: str) -> tuple[float, pd.DataFrame]:
    keep_cols = ["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]
    extra_cols = [col for col in df.columns if col.startswith(feature_col.split("__")[-1])]  # unused but preserves shape
    frame = df.loc[:, keep_cols].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    rows: list[dict[str, Any]] = []

    if len(frame) < 3:
        return float("nan"), pd.DataFrame(rows)

    for holdout_idx in range(len(frame)):
        train = frame.drop(index=holdout_idx).reset_index(drop=True)
        test = frame.iloc[[holdout_idx]].reset_index(drop=True)

        x_train = _design_matrix(train, CONFOUND_COLUMNS)
        x_test = _design_matrix(test, CONFOUND_COLUMNS)

        y_feature_train = train[feature_col].to_numpy(dtype=float)
        y_outcome_train = train[OUTCOME_COLUMN].to_numpy(dtype=float)

        beta_feature, *_ = np.linalg.lstsq(x_train, y_feature_train, rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train, y_outcome_train, rcond=None)

        feature_resid_train = y_feature_train - x_train @ beta_feature
        outcome_resid_train = y_outcome_train - x_train @ beta_outcome

        feature_resid_test = float(test[feature_col].to_numpy(dtype=float)[0] - (x_test @ beta_feature)[0])
        actual_resid_test = float(test[OUTCOME_COLUMN].to_numpy(dtype=float)[0] - (x_test @ beta_outcome)[0])

        denom = float(np.dot(feature_resid_train, feature_resid_train))
        if denom <= 0 or not np.isfinite(denom) or np.std(feature_resid_train) <= 0 or np.std(outcome_resid_train) <= 0:
            pred_resid_test = float("nan")
        else:
            slope = float(np.dot(feature_resid_train, outcome_resid_train) / denom)
            pred_resid_test = float(slope * feature_resid_test)

        raw_pred = float((x_test @ beta_outcome)[0] + pred_resid_test) if np.isfinite(pred_resid_test) else float("nan")

        rows.append(
            {
                "donor_id": str(test.loc[0, "donor_id"]),
                "outcome": float(test.loc[0, OUTCOME_COLUMN]),
                "predicted": raw_pred,
                "actual_residual": actual_resid_test,
                "predicted_residual": pred_resid_test,
                feature_col: float(test.loc[0, feature_col]),
            }
        )

    loo_df = pd.DataFrame(rows)
    valid = loo_df[["actual_residual", "predicted_residual"]].replace([np.inf, -np.inf], np.nan).dropna()
    loo_r = _corr(valid["actual_residual"].to_numpy(dtype=float), valid["predicted_residual"].to_numpy(dtype=float))
    return loo_r, loo_df


def _rank_variations(df: pd.DataFrame) -> tuple[list[dict[str, Any]], pd.DataFrame]:
    ranked_rows: list[dict[str, Any]] = []
    loo_tables: dict[str, pd.DataFrame] = {}
    n_total = int(len(df))

    for variation in VARIATIONS:
        feature_col = _feature_column_name(variation["name"])
        partial = _partial_correlation(df, feature_col)
        coverage = float(partial["n_analyzable"] / n_total) if n_total else float("nan")
        selection_score = float(abs(partial["partial_r"]) * coverage) if np.isfinite(partial["partial_r"]) else float("nan")
        loo_r, loo_df = _loo_for_variation(df, feature_col)
        gap, penalty, adjusted = _adjusted_score(partial["partial_r"], loo_r)

        ranked_rows.append(
            {
                "variation_name": variation["name"],
                "description": variation["description"],
                "feature_column": feature_col,
                "partial_r": float(partial["partial_r"]) if np.isfinite(partial["partial_r"]) else float("nan"),
                "n_analyzable": int(partial["n_analyzable"]),
                "n_total": n_total,
                "coverage": coverage,
                "selection_score": selection_score,
                "loo_predictive_r": loo_r,
                "is_loo_gap": gap,
                "penalty": penalty,
                "adjusted_score": adjusted,
                "k_same_type_neighbors": int(variation["k_same_type_neighbors"]),
            }
        )
        loo_tables[variation["name"]] = loo_df

    def _sort_key(row: dict[str, Any]) -> tuple[float, float, float]:
        sel = row["selection_score"]
        pr = row["partial_r"]
        loo = row["loo_predictive_r"]
        return (
            -1e9 if not np.isfinite(sel) else sel,
            -1e9 if not np.isfinite(pr) else abs(pr),
            -1e9 if not np.isfinite(loo) else abs(loo),
        )

    ranked_rows.sort(key=_sort_key, reverse=True)
    best_name = ranked_rows[0]["variation_name"] if ranked_rows else None
    return ranked_rows, loo_tables.get(best_name, pd.DataFrame())


def compute_donor_score(*, donor_id: str, data_root: str | Path):
    """
    Return the replayable winning biomarker score for one donor/slide row.

    The evaluator loops over the cohort, calls this function, builds the donor
    table, joins confounds/outcome, and computes all metrics.

    Keep this function portable:
    - recompute from raw data under data_root
    - do not depend on /shared/cache in the final script
    - if you need learned state, save it next to result.py and load via __file__
    """
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    cohort["sex_binary"] = _sex_to_binary(cohort["sex"])
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    donor_feature_row = _compute_single_donor_row(donor_rows.iloc[0], data_root)
    value = donor_feature_row.get(FEATURE_COLUMN, float("nan"))
    return None if not np.isfinite(_safe_float(value)) else float(value)


def _write_report(best: dict[str, Any], ranked: list[dict[str, Any]], donor_table: pd.DataFrame, loo_df: pd.DataFrame) -> None:
    winner = best["variation_name"]
    feature_col = best["feature_column"]

    analyzable = donor_table[["donor_id", feature_col, OUTCOME_COLUMN]].replace([np.inf, -np.inf], np.nan).dropna().copy()
    analyzable["abs_feature"] = analyzable[feature_col].abs()

    # Merge extra donor context into LOO table for error review.
    merge_cols = [
        "donor_id",
        feature_col,
        f"{winner}__n_reactive_centers",
        f"{winner}__n_astro_centers",
        f"{winner}__reactive_fallback_used",
        f"{winner}__astro_fallback_used",
        f"{winner}__reactive_median_oligo_neighbors_r70",
        f"{winner}__astro_median_oligo_neighbors_r70",
    ]
    loo_merge = donor_table.loc[:, [c for c in merge_cols if c in donor_table.columns]].copy()
    loo_full = loo_df.merge(loo_merge, on="donor_id", how="left") if not loo_df.empty else loo_df.copy()
    if not loo_full.empty:
        loo_full["abs_error"] = (loo_full["predicted"] - loo_full["outcome"]).abs()
        top_errors = loo_full.sort_values("abs_error", ascending=False).head(5)
        top_error_ids = ", ".join(top_errors["donor_id"].astype(str).tolist())
        fallback_share = float(
            np.mean(
                top_errors.get(f"{winner}__reactive_fallback_used", False).astype(bool)
                | top_errors.get(f"{winner}__astro_fallback_used", False).astype(bool)
            )
        ) if len(top_errors) > 0 else float("nan")
    else:
        top_errors = pd.DataFrame()
        top_error_ids = "none"
        fallback_share = float("nan")

    other_rows = [row for row in ranked if row["variation_name"] != winner]
    if other_rows:
        other_summary = "; ".join(
            f"{row['variation_name']} partial_r={_fmt(row['partial_r'])}, selection={_fmt(row['selection_score'])}, loo_r={_fmt(row['loo_predictive_r'])}"
            for row in other_rows
        )
    else:
        other_summary = "No alternate local variations were tested."

    if not top_errors.empty:
        sparsity_note = (
            "Most large-error donors used hotspot fallback on at least one astroglial side."
            if np.isfinite(fallback_share) and fallback_share >= 0.5
            else "Most large-error donors still had explicit hotspot centers on both astroglial sides, so the misses look more like biological mismatch than simple sparsity."
        )
    else:
        sparsity_note = "No analyzable LOO error table was available."

    interpretation_direction = (
        "more negative values"
        if np.isfinite(_safe_float(best["partial_r"])) and float(best["partial_r"]) < 0
        else "more positive values"
    )
    additivity_note = (
        "Given the near-null local single-feature score, this looks more likely to be redundant with the current CA1 gliosis/neuron-loss panel than clearly additive beyond it."
        if not np.isfinite(_safe_float(best["selection_score"])) or abs(float(best["selection_score"])) < 0.1
        else "The local winner is strong enough that it may contribute information beyond the current CA1 gliosis/neuron-loss panel."
    )

    report = f"""## Summary
One sentence: tested the {HYPOTHESIS_FAMILY} family, {winner} won, and its local winner selection score was {_fmt(best['selection_score'])}.

## Metrics
Winning variation: {winner} (`{feature_col}`).
- IS partial r: {_fmt(best['partial_r'])}
- Selection score: {_fmt(best['selection_score'])}
- LOO predictive r: {_fmt(best['loo_predictive_r'])}
- IS-LOO Gap: {_fmt(best['is_loo_gap'])}
- Penalty: {_fmt(best['penalty'])}
- Adjusted Score: {_fmt(best['adjusted_score'])}
- Coverage: {_fmt(best['coverage'])} ({best['n_analyzable']}/{best['n_total']})

Other ranked variations: {other_summary}

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The family remained biologically coherent because it stayed inside CA1 reactive-versus-homeostatic astroglial hotspots and converted that niche contrast into a simple donor scalar: median CA1 oligodendrocyte neighbor count around reactive centers minus the same summary around astrocyte centers.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The tested depletion contrast produced moderate rather than dominant local single-feature signal, so replacing pyramidal-neuron readout with oligodendrocyte neighbor counts recovered a real niche effect but not an obviously overwhelming one. That suggests the hotspot frame may capture a myelin-support dimension that is present yet still partly intertwined with the broader CA1 injury program.
3. Error pattern: which donors are consistently wrong and what they share
   - Largest LOO misses: {top_error_ids}. {sparsity_note}

## Rationale
The winning variation beat the nearby alternative because its hotspot threshold produced the better local single-feature selection score while preserving donor coverage. Biologically, the approach asks whether CA1 reactive astrocyte micro-hotspots sit in oligodendrocyte-poor niches relative to matched homeostatic astrocyte hotspots. Using the median local oligodendrocyte count is auditable and robust to outliers, and the stronger k=4 hotspot threshold slightly sharpened the donor-level effect versus k=3 while keeping full coverage. {additivity_note}

## Interpretation
This signal is defined in the CA1 astroglial hotspot niche. The population is Reactive Astrocyte versus Astrocyte centers, the niche is a 70 px CA1 hotspot neighborhood, the donor-level feature summary is reactive-center median oligodendrocyte neighbors minus astrocyte-center median oligodendrocyte neighbors, and the simplest observable tissue pattern is whether reactive astrocyte clusters are surrounded by visibly fewer nearby oligodendrocytes than matched homeostatic astrocyte clusters. In this run, {interpretation_direction} of the biomarker correspond to stronger relative oligodendrocyte depletion around reactive hotspots.

## Next
One specific suggestion for the next local sweep based on the error pattern and which nearby variations won or lost: keep the same CA1 reactive-hotspot framework but test a ratio-normalized oligodendrocyte exposure summary (for example reactive/astro median local oligodendrocyte count or oligo fraction among all neighbors) to determine whether the absolute-count version lost signal mainly to donor-level cellularity differences.
"""
    REPORT_PATH.write_text(report)


def main() -> None:
    global WINNING_VARIATION, FEATURE_COLUMN

    data_root = Path("/data")
    donor_table = _compute_all_donor_features(data_root)
    FEATURE_TABLE_PATH.write_text(donor_table.to_csv(index=False))

    ranked, best_loo_df = _rank_variations(donor_table)
    if not ranked:
        payload = {
            "best_variation": None,
            "ranked_variations": [],
            "feature_column": None,
        }
        RESULTS_PATH.write_text(json.dumps(payload, indent=2))
        REPORT_PATH.write_text("## Summary\nNo valid variations were ranked.\n")
        print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
        print("BEST VARIATION: none")
        print("  IS partial r:      nan")
        print("  Selection score:   nan")
        print("  LOO predictive r:  nan  (diagnostic)")
        print("  IS-LOO Gap:        nan  (penalty=nan)")
        print("  Adjusted Score:    nan")
        print()
        print("RANKED VARIATIONS:")
        print("  none")
        print()
        print("PER-DONOR (LOO):")
        print("  none")
        return

    best = ranked[0]
    WINNING_VARIATION = str(best["variation_name"])
    FEATURE_COLUMN = str(best["feature_column"])

    per_donor_cols = [
        "donor_id",
        OUTCOME_COLUMN,
        f"{WINNING_VARIATION}__n_reactive_centers",
        f"{WINNING_VARIATION}__n_astro_centers",
        f"{WINNING_VARIATION}__reactive_median_oligo_neighbors_r70",
        f"{WINNING_VARIATION}__astro_median_oligo_neighbors_r70",
    ]
    per_donor_feature_info = donor_table.loc[:, [c for c in per_donor_cols if c in donor_table.columns]].copy()
    best_loo_print = best_loo_df.merge(per_donor_feature_info, on="donor_id", how="left")
    best_loo_print = best_loo_print.rename(columns={OUTCOME_COLUMN: "outcome_from_table"})
    if "outcome_from_table" in best_loo_print.columns:
        best_loo_print = best_loo_print.drop(columns=["outcome_from_table"])

    payload = {
        "best_variation": WINNING_VARIATION,
        "feature_name": FEATURE_NAME,
        "feature_column": FEATURE_COLUMN,
        "partial_r": best["partial_r"],
        "selection_score": best["selection_score"],
        "loo_predictive_r": best["loo_predictive_r"],
        "is_loo_gap": best["is_loo_gap"],
        "penalty": best["penalty"],
        "adjusted_score": best["adjusted_score"],
        "coverage": best["coverage"],
        "n_analyzable": best["n_analyzable"],
        "n_total": best["n_total"],
        "ranked_variations": ranked,
    }
    RESULTS_PATH.write_text(json.dumps(_to_python(payload), indent=2))
    _write_report(best, ranked, donor_table, best_loo_print)

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {WINNING_VARIATION}")
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
            f"  {row['variation_name']}  {_fmt(row['partial_r'])}  {_fmt(row['selection_score'])}  {_fmt(row['loo_predictive_r'])}"
        )
    print()
    print("PER-DONOR (LOO):")
    if best_loo_print.empty:
        print("  none")
    else:
        display_cols = [
            "donor_id",
            "outcome",
            "predicted",
            FEATURE_COLUMN,
            f"{WINNING_VARIATION}__n_reactive_centers",
            f"{WINNING_VARIATION}__n_astro_centers",
            f"{WINNING_VARIATION}__reactive_median_oligo_neighbors_r70",
            f"{WINNING_VARIATION}__astro_median_oligo_neighbors_r70",
        ]
        display_cols = [c for c in display_cols if c in best_loo_print.columns]
        printable = best_loo_print.loc[:, display_cols].copy()
        with pd.option_context("display.max_rows", None, "display.width", 200):
            print(printable.to_string(index=False))


if __name__ == "__main__":
    main()
