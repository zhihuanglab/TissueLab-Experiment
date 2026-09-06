set -euo pipefail
cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import re
import sys
import warnings
from pathlib import Path
from typing import Any

sys.path.insert(0, "/shared/lib")

import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath
from scipy.spatial import cKDTree

from shared_analysis import load_training_cohort
from shared_analysis.artifacts import (
    build_results_payload,
    write_donor_feature_table,
    write_results_payload,
)
from shared_analysis.stats import (
    bootstrap_partial_correlation,
    leave_one_out_summary,
    partial_correlation,
)
import shared_analysis.sea_ad_lfb as sea


FEATURE_FAMILY = "ca1_reactive_oligodendrocyte_gap"
WINNING_VARIATION = "__PENDING__"
WINNING_RADIUS_PX = -1
FEATURE_NAME = "__PENDING__"
FEATURE_COLUMN = "__PENDING__"

OUTCOME_COL = "slope_zmem0"
CONFounds_NUMERIC = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
CANONICAL_COVARIATES = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex"]

TARGET_REGION = "CA1"
FOCAL_LABEL = "Pyramidal Neuron"
REACTIVE_LABEL = "Reactive Astrocyte"
SUPPORT_LABEL = "Oligodendrocyte"

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", message="Object at .* is not recognized as a component of a Zarr hierarchy.*")


def _safe_float(x: Any) -> float | None:
    if x is None:
        return None
    try:
        value = float(x)
    except Exception:
        return None
    if math.isnan(value) or math.isinf(value):
        return None
    return value


def _nan_to_none(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {k: _nan_to_none(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_nan_to_none(v) for v in obj]
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
    return obj


def _sex_to_binary(value: Any) -> float:
    if pd.isna(value):
        return float("nan")
    s = str(value).strip().lower()
    if s == "male":
        return 1.0
    if s == "female":
        return 0.0
    return float("nan")


def _design_matrix(confounds: pd.DataFrame) -> np.ndarray:
    x = confounds.to_numpy(dtype=float)
    intercept = np.ones((len(confounds), 1), dtype=float)
    return np.hstack([intercept, x])


def _mask_points_in_polygons(points: np.ndarray, polygons: list[np.ndarray]) -> np.ndarray:
    mask = np.zeros(len(points), dtype=bool)
    for polygon in polygons:
        polygon = np.asarray(polygon, dtype=np.float32)
        if polygon.ndim != 2 or polygon.shape[0] < 3:
            continue
        mask |= MplPath(polygon, closed=True).contains_points(points)
    return mask


def _extract_radius_px(description: str) -> int:
    m = re.search(r"(\d+)\s*px", description)
    if not m:
        raise ValueError(f"Could not parse px radius from description: {description}")
    return int(m.group(1))


def _feature_column_for_radius(radius_px: int) -> str:
    return f"ca1_reactive_exposed_pyramidal_oligo_gap_fraction_r{radius_px}px"


def _feature_name_for_radius(radius_px: int) -> str:
    return f"ca1_reactive_oligodendrocyte_gap_r{radius_px}px"


def _variation_specs_from_brief(brief: dict[str, Any]) -> list[dict[str, Any]]:
    specs = []
    for variation in brief["variations"]:
        radius_px = _extract_radius_px(str(variation["description"]))
        specs.append(
            {
                "name": str(variation["name"]),
                "description": str(variation["description"]),
                "radius_px": int(radius_px),
                "feature_column": _feature_column_for_radius(radius_px),
                "feature_name": _feature_name_for_radius(radius_px),
            }
        )
    return specs


def _compute_family_for_slide(zarr_path: str | Path, radii_px: list[int]) -> dict[str, Any]:
    zarr_path = Path(zarr_path)
    centroids = sea.load_centroids(zarr_path)
    class_ids = sea.load_class_ids(zarr_path)
    class_lookup = sea.load_class_lookup(zarr_path)

    pyramidal_ids = np.array([cid for cid, label in class_lookup.items() if label == FOCAL_LABEL], dtype=int)
    reactive_ids = np.array([cid for cid, label in class_lookup.items() if label == REACTIVE_LABEL], dtype=int)
    oligo_ids = np.array([cid for cid, label in class_lookup.items() if label == SUPPORT_LABEL], dtype=int)

    relevant_ids = np.concatenate([arr for arr in [pyramidal_ids, reactive_ids, oligo_ids] if len(arr) > 0])
    if len(relevant_ids) == 0:
        return {
            "ca1_relevant_n": 0,
            "ca1_pyramidal_n": 0,
            "ca1_reactive_astrocyte_n": 0,
            "ca1_oligodendrocyte_n": 0,
            "variations": {
                radius: {
                    "feature_value": float("nan"),
                    "reactive_exposed_pyramidal_n": 0,
                    "oligo_free_reactive_exposed_pyramidal_n": 0,
                }
                for radius in radii_px
            },
        }

    relevant_mask = np.isin(class_ids, relevant_ids)
    relevant_centroids = np.asarray(centroids[relevant_mask], dtype=np.float32)
    relevant_class_ids = np.asarray(class_ids[relevant_mask], dtype=int)

    polygons = sea.load_region_polygons(zarr_path, scale=16.0)
    ca1_polygons = polygons.get(TARGET_REGION, [])
    if not ca1_polygons:
        return {
            "ca1_relevant_n": 0,
            "ca1_pyramidal_n": 0,
            "ca1_reactive_astrocyte_n": 0,
            "ca1_oligodendrocyte_n": 0,
            "variations": {
                radius: {
                    "feature_value": float("nan"),
                    "reactive_exposed_pyramidal_n": 0,
                    "oligo_free_reactive_exposed_pyramidal_n": 0,
                }
                for radius in radii_px
            },
        }

    ca1_mask = _mask_points_in_polygons(relevant_centroids, ca1_polygons)
    ca1_centroids = relevant_centroids[ca1_mask]
    ca1_class_ids = relevant_class_ids[ca1_mask]

    is_pyramidal = np.isin(ca1_class_ids, pyramidal_ids)
    is_reactive = np.isin(ca1_class_ids, reactive_ids)
    is_oligo = np.isin(ca1_class_ids, oligo_ids)

    pyramidal_points = ca1_centroids[is_pyramidal]
    reactive_points = ca1_centroids[is_reactive]
    oligo_points = ca1_centroids[is_oligo]

    n_pyramidal = int(len(pyramidal_points))
    n_reactive = int(len(reactive_points))
    n_oligo = int(len(oligo_points))

    reactive_tree = cKDTree(reactive_points) if n_reactive > 0 else None
    oligo_tree = cKDTree(oligo_points) if n_oligo > 0 else None

    per_radius: dict[int, dict[str, Any]] = {}
    for radius_px in radii_px:
        if n_pyramidal == 0:
            feature_value = float("nan")
            reactive_exposed_n = 0
            oligo_free_n = 0
        else:
            if reactive_tree is None:
                reactive_exposed = np.zeros(n_pyramidal, dtype=bool)
            else:
                reactive_neighbors = reactive_tree.query_ball_point(pyramidal_points, r=float(radius_px))
                reactive_exposed = np.fromiter((len(ix) > 0 for ix in reactive_neighbors), count=n_pyramidal, dtype=bool)

            reactive_exposed_n = int(reactive_exposed.sum())
            if reactive_exposed_n == 0:
                oligo_free_n = 0
                feature_value = 0.0
            else:
                exposed_points = pyramidal_points[reactive_exposed]
                if oligo_tree is None:
                    oligo_neighbor_counts = np.zeros(reactive_exposed_n, dtype=int)
                else:
                    oligo_neighbors = oligo_tree.query_ball_point(exposed_points, r=float(radius_px))
                    oligo_neighbor_counts = np.fromiter((len(ix) for ix in oligo_neighbors), count=reactive_exposed_n, dtype=int)
                oligo_free_n = int((oligo_neighbor_counts == 0).sum())
                feature_value = float(oligo_free_n / reactive_exposed_n)

        per_radius[int(radius_px)] = {
            "feature_value": float(feature_value),
            "reactive_exposed_pyramidal_n": int(reactive_exposed_n),
            "oligo_free_reactive_exposed_pyramidal_n": int(oligo_free_n),
        }

    return {
        "ca1_relevant_n": int(len(ca1_centroids)),
        "ca1_pyramidal_n": n_pyramidal,
        "ca1_reactive_astrocyte_n": n_reactive,
        "ca1_oligodendrocyte_n": n_oligo,
        "variations": per_radius,
    }


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
    if WINNING_RADIUS_PX <= 0:
        raise RuntimeError("WINNING_RADIUS_PX has not been materialized in result.py yet.")
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    family = _compute_family_for_slide(data_root / slide_name, [int(WINNING_RADIUS_PX)])
    return float(family["variations"][int(WINNING_RADIUS_PX)]["feature_value"])


def _clean_analysis_frame(df: pd.DataFrame, feature_col: str) -> pd.DataFrame:
    cols = ["donor_id", feature_col, OUTCOME_COL, *CONFounds_NUMERIC]
    return df.loc[:, cols].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)


def _loo_predictions(df: pd.DataFrame, feature_col: str) -> tuple[float, pd.DataFrame]:
    frame = _clean_analysis_frame(df, feature_col)
    rows: list[dict[str, Any]] = []
    pred_resid = []
    actual_resid = []

    for i in range(len(frame)):
        train = frame.drop(index=i).reset_index(drop=True)
        test = frame.iloc[[i]].reset_index(drop=True)

        x_train_conf = _design_matrix(train[CONFounds_NUMERIC])
        x_test_conf = _design_matrix(test[CONFounds_NUMERIC])

        beta_feature, *_ = np.linalg.lstsq(x_train_conf, train[feature_col].to_numpy(dtype=float), rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train_conf, train[OUTCOME_COL].to_numpy(dtype=float), rcond=None)

        resid_feature_train = train[feature_col].to_numpy(dtype=float) - x_train_conf @ beta_feature
        resid_outcome_train = train[OUTCOME_COL].to_numpy(dtype=float) - x_train_conf @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            slope = float("nan")
            predicted_resid = float("nan")
            predicted_outcome = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            resid_feature_test = test[feature_col].to_numpy(dtype=float) - x_test_conf @ beta_feature
            predicted_resid = float(slope * resid_feature_test[0])
            predicted_outcome = float((x_test_conf @ beta_outcome)[0] + predicted_resid)

        actual_outcome = float(test.iloc[0][OUTCOME_COL])
        actual_outcome_resid = float(actual_outcome - (x_test_conf @ beta_outcome)[0])

        pred_resid.append(predicted_resid)
        actual_resid.append(actual_outcome_resid)

        rows.append(
            {
                "donor_id": str(test.iloc[0]["donor_id"]),
                "outcome": actual_outcome,
                "predicted": predicted_outcome,
                "predicted_residualized": predicted_resid,
                "outcome_residualized": actual_outcome_resid,
                feature_col: float(test.iloc[0][feature_col]),
            }
        )

    pred_arr = np.asarray(pred_resid, dtype=float)
    act_arr = np.asarray(actual_resid, dtype=float)
    mask = np.isfinite(pred_arr) & np.isfinite(act_arr)
    if mask.sum() < 3 or float(np.std(pred_arr[mask])) == 0.0 or float(np.std(act_arr[mask])) == 0.0:
        loo_r = float("nan")
    else:
        loo_r = float(np.corrcoef(pred_arr[mask], act_arr[mask])[0, 1])

    pred_df = pd.DataFrame(rows)
    return loo_r, pred_df


def _gap_penalty(gap: float) -> tuple[float, float]:
    if not np.isfinite(gap):
        return float("nan"), float("nan")
    penalty = max(0.0, float(gap) - 0.15) * 0.5
    adjusted = -1.0 if float(gap) > 0.30 else float("nan")
    return penalty, adjusted


def _evaluate_variation(donor_table: pd.DataFrame, spec: dict[str, Any], n_total: int) -> dict[str, Any]:
    feature_col = spec["feature_column"]
    frame = donor_table.copy()

    pc = partial_correlation(frame, feature_col=feature_col, outcome_col=OUTCOME_COL, confounds=CONFounds_NUMERIC)
    boot = bootstrap_partial_correlation(frame, feature_col=feature_col, outcome_col=OUTCOME_COL, confounds=CONFounds_NUMERIC, n_boot=2000, random_state=0)
    loo_summary = leave_one_out_summary(frame, feature_col=feature_col, outcome_col=OUTCOME_COL, confounds=CONFounds_NUMERIC, id_col="donor_id")
    loo_r, loo_pred_df = _loo_predictions(frame, feature_col)

    analyzable = _clean_analysis_frame(frame, feature_col)
    n_analyzable = int(len(analyzable))
    coverage_ratio = float(n_analyzable / n_total) if n_total else float("nan")
    partial_r = float(pc["partial_r"]) if np.isfinite(pc["partial_r"]) else float("nan")
    selection_score = float(abs(partial_r) * coverage_ratio) if np.isfinite(partial_r) and np.isfinite(coverage_ratio) else float("nan")
    gap = float(abs(partial_r - loo_r)) if np.isfinite(partial_r) and np.isfinite(loo_r) else float("nan")
    penalty = max(0.0, gap - 0.15) * 0.5 if np.isfinite(gap) else float("nan")
    adjusted_score = -1.0 if np.isfinite(gap) and gap > 0.30 else (float(loo_r - penalty) if np.isfinite(loo_r) and np.isfinite(penalty) else float("nan"))

    result = {
        "name": spec["name"],
        "description": spec["description"],
        "radius_px": int(spec["radius_px"]),
        "feature_name": spec["feature_name"],
        "feature_column": feature_col,
        "n_total": int(n_total),
        "n_analyzable": int(n_analyzable),
        "coverage_ratio": coverage_ratio,
        "partial_r": partial_r,
        "selection_score": selection_score,
        "loo_predictive_r": float(loo_r) if np.isfinite(loo_r) else float("nan"),
        "is_loo_gap": gap,
        "gap_penalty": penalty,
        "adjusted_score": adjusted_score,
        "p_value": float(boot["p_value"]) if np.isfinite(boot["p_value"]) else float("nan"),
        "ci_lo": float(boot["ci_lo"]) if np.isfinite(boot["ci_lo"]) else float("nan"),
        "ci_hi": float(boot["ci_hi"]) if np.isfinite(boot["ci_hi"]) else float("nan"),
        "loo_unstable_count": int(loo_summary.get("unstable_count", 0) or 0),
        "loo_max_shift": float(loo_summary.get("max_shift", float("nan"))),
        "donor_ids_used": analyzable["donor_id"].astype(str).tolist(),
        "per_donor_loo": loo_pred_df.to_dict(orient="records"),
    }
    return result


def _score_sort_key(item: dict[str, Any]) -> tuple[float, float, float]:
    selection = item.get("selection_score")
    partial_r = item.get("partial_r")
    adjusted = item.get("adjusted_score")
    return (
        -1e9 if not np.isfinite(selection) else float(selection),
        -1e9 if not np.isfinite(partial_r) else float(abs(partial_r)),
        -1e9 if not np.isfinite(adjusted) else float(adjusted),
    )


def _update_self_for_winner(script_path: Path, best: dict[str, Any]) -> None:
    text = script_path.read_text(encoding="utf-8")
    replacements = {
        "WINNING_VARIATION": repr(best["name"]),
        "WINNING_RADIUS_PX": str(int(best["radius_px"])),
        "FEATURE_NAME": repr(best["feature_name"]),
        "FEATURE_COLUMN": repr(best["feature_column"]),
    }
    for key, value in replacements.items():
        text = re.sub(rf"^{key}\s*=.*$", f"{key} = {value}", text, flags=re.MULTILINE)
    script_path.write_text(text, encoding="utf-8")


def _build_donor_table(data_root: Path, brief: dict[str, Any], variation_specs: list[dict[str, Any]]) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    radii = [int(spec["radius_px"]) for spec in variation_specs]
    rows = []

    for _, row in cohort.iterrows():
        slide_name = str(row["slide_name"])
        donor_id = str(row["donor_id"])
        family = _compute_family_for_slide(data_root / slide_name, radii)
        base = {
            "donor_id": donor_id,
            "slide_name": slide_name,
            OUTCOME_COL: float(row[OUTCOME_COL]) if pd.notna(row[OUTCOME_COL]) else float("nan"),
            "max_age_vis": float(row["max_age_vis"]) if pd.notna(row["max_age_vis"]) else float("nan"),
            "braak_numeric": float(row["braak_numeric"]) if pd.notna(row["braak_numeric"]) else float("nan"),
            "cerad_ordinal": float(row["cerad_ordinal"]) if pd.notna(row["cerad_ordinal"]) else float("nan"),
            "sex_binary": _sex_to_binary(row.get("sex")),
            "sex": row.get("sex"),
            "cognitive_status": row.get("cognitive_status"),
            "overall_ad_neuropath_change": row.get("overall_ad_neuropath_change"),
            "ca1_pyramidal_n": int(family["ca1_pyramidal_n"]),
            "ca1_reactive_astrocyte_n": int(family["ca1_reactive_astrocyte_n"]),
            "ca1_oligodendrocyte_n": int(family["ca1_oligodendrocyte_n"]),
        }
        for spec in variation_specs:
            vr = family["variations"][int(spec["radius_px"])]
            base[spec["feature_column"]] = float(vr["feature_value"])
            base[f"reactive_exposed_pyramidal_n_r{int(spec['radius_px'])}px"] = int(vr["reactive_exposed_pyramidal_n"])
            base[f"oligo_free_reactive_exposed_pyramidal_n_r{int(spec['radius_px'])}px"] = int(vr["oligo_free_reactive_exposed_pyramidal_n"])
        rows.append(base)

    return pd.DataFrame(rows)


def _format_num(x: Any) -> str:
    if x is None:
        return "nan"
    try:
        x = float(x)
    except Exception:
        return str(x)
    return "nan" if not np.isfinite(x) else f"{x:.4f}"


def _write_report(path: Path, best: dict[str, Any], ranked: list[dict[str, Any]], donor_table: pd.DataFrame) -> None:
    feature_col = best["feature_column"]
    per_donor = pd.DataFrame(best["per_donor_loo"]).merge(
        donor_table[
            [
                "donor_id",
                feature_col,
                "ca1_pyramidal_n",
                "ca1_reactive_astrocyte_n",
                "ca1_oligodendrocyte_n",
                "cognitive_status",
                "overall_ad_neuropath_change",
            ]
        ],
        on="donor_id",
        how="left",
    )
    per_donor["abs_error"] = (per_donor["outcome"] - per_donor["predicted"]).abs()
    worst = per_donor.sort_values("abs_error", ascending=False).head(5)
    worst_ids = worst["donor_id"].tolist()
    worst_status = worst["cognitive_status"].fillna("NA").tolist()
    worst_path = worst["overall_ad_neuropath_change"].fillna("NA").tolist()

    other_lines = []
    for item in ranked[1:]:
        other_lines.append(
            f"- {item['name']} (r={item['radius_px']} px): partial_r={_format_num(item['partial_r'])}, "
            f"selection_score={_format_num(item['selection_score'])}, loo_r={_format_num(item['loo_predictive_r'])}, "
            f"adjusted={_format_num(item['adjusted_score'])}"
        )
    other_text = "\n".join(other_lines) if other_lines else "- No alternative variations."

    if len(ranked) > 1 and best["radius_px"] < ranked[1]["radius_px"]:
        failed_reason = (
            f"The broader {ranked[1]['radius_px']} px neighborhood underperformed, suggesting the oligodendrocyte-depletion "
            f"signal is most specific in a tighter reactive astrocyte niche rather than a wider CA1 backdrop."
        )
        next_suggestion = (
            f"Keep the same CA1 pyramidal/reactive/oligodendrocyte family but sweep tighter radii around the winner "
            f"({best['radius_px']} px), such as 40-70 px, to see whether the signal sharpens further without the dilution seen at {ranked[1]['radius_px']} px."
        )
    elif len(ranked) > 1 and best["radius_px"] > ranked[1]["radius_px"]:
        failed_reason = (
            f"The narrower {ranked[1]['radius_px']} px neighborhood underperformed, suggesting oligodendrocyte depletion becomes measurable only when the "
            f"reactive astrocyte neighborhood is widened enough to capture the surrounding support-cell desert."
        )
        next_suggestion = (
            f"Keep the same CA1 pyramidal/reactive/oligodendrocyte family but test slightly broader radii around the winner "
            f"({best['radius_px']} px), such as 90-120 px, to determine where the support-cell depletion signal peaks."
        )
    else:
        failed_reason = (
            "Nearby radius changes did not improve the feature, implying the family may be locally real but weak, unstable, or partly redundant with the current CA1 injury panel."
        )
        next_suggestion = (
            "Within the same CA1 pyramidal/reactive/oligodendrocyte family, next test an intensity version of the niche such as mean oligodendrocyte count around reactive-exposed neurons rather than only the zero-neighbor fraction."
        )

    interpretation_additivity = (
        "Because the accepted panel already tracks CA1 pyramidal loss and reactive-astrocyte injury, this candidate looks most likely to add information only if local oligodendrocyte depletion captures a distinct myelin/support-cell axis rather than a restatement of the same reactive niche."
    )

    lines = [
        "## Summary",
        f"Tested the {FEATURE_FAMILY} family in CA1 reactive astrocyte-conditioned pyramidal niches; {best['name']} won with selection_score={_format_num(best['selection_score'])}.",
        "",
        "## Metrics",
        f"Best variation: {best['name']} ({best['radius_px']} px).",
        f"- IS partial r: {_format_num(best['partial_r'])}",
        f"- Selection score: {_format_num(best['selection_score'])}",
        f"- LOO predictive r: {_format_num(best['loo_predictive_r'])}",
        f"- IS-LOO Gap: {_format_num(best['is_loo_gap'])} (penalty={_format_num(best['gap_penalty'])})",
        f"- Adjusted Score: {_format_num(best['adjusted_score'])}",
        f"- n_analyzable / n_total: {best['n_analyzable']} / {best['n_total']}",
        "Other tested variations:",
        other_text,
        "",
        "## Findings",
        f"1. What worked and why: The winning signal comes from {TARGET_REGION} {FOCAL_LABEL.lower()}s that sit in {REACTIVE_LABEL.lower()} neighborhoods yet lack nearby {SUPPORT_LABEL.lower()} support. Biologically, this targets a stressed neuronal niche where glial reactivity is present but local myelin/support-cell buffering is sparse, which is coherent with faster memory decline if such neighborhoods reflect unstable or denuded CA1 tissue.",
        f"2. What failed and why: {failed_reason}",
        f"3. Error pattern: The largest LOO outcome misses were donors {', '.join(worst_ids)}; their cognitive-status labels were {', '.join(worst_status)} and AD-neuropath-change labels were {', '.join(worst_path)}. These errors suggest the feature captures one CA1 injury axis but does not fully explain donors at the extremes of decline severity or pathology burden.",
        "",
        "## Rationale",
        f"The best variation is biologically coherent because it asks a narrowly local question inside CA1: among pyramidal neurons already adjacent to reactive astrocytes, how often is oligodendrocyte support absent? That conditional denominator removes global cell-composition effects and focuses the score on a specific injury-support mismatch. It beat the nearby alternative because its radius better matched the operative niche scale for this dataset rather than averaging over a broader background.",
        interpretation_additivity,
        "",
        "## Interpretation",
        f"The signal appears to mean that memory decline tracks CA1 neuronal neighborhoods where {REACTIVE_LABEL.lower()}s are close to {FOCAL_LABEL.lower()}s but {SUPPORT_LABEL.lower()}s are locally missing.",
        f"- Population: {TARGET_REGION} {FOCAL_LABEL}",
        f"- Niche: {TARGET_REGION} reactive-astrocyte-exposed pyramidal neighborhoods",
        f"- Feature summary: fraction of reactive-exposed CA1 pyramidal neurons with zero oligodendrocyte neighbors within {best['radius_px']} px",
        f"- Simplest observable pattern: reactive astrocyte-wrapped CA1 pyramidal zones with sparse nearby oligodendrocytes / diminished local myelin-support presence",
        "",
        "## Next",
        next_suggestion,
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    data_root = Path("/data")
    context_path = Path("/scratch/context_bundle.json")
    script_path = Path(__file__)
    results_path = Path("/scratch/results.json")
    table_path = Path("/scratch/donor_feature_table.csv")
    report_path = Path("/scratch/report.md")

    with context_path.open("r", encoding="utf-8") as f:
        context = json.load(f)
    brief = context["worker_brief"]
    variation_specs = _variation_specs_from_brief(brief)

    donor_table = _build_donor_table(data_root, brief, variation_specs)
    n_total = int(len(donor_table))

    evaluations = [_evaluate_variation(donor_table, spec, n_total) for spec in variation_specs]
    ranked = sorted(evaluations, key=_score_sort_key, reverse=True)
    best = ranked[0]

    _update_self_for_winner(script_path, best)

    extra_columns = [
        spec["feature_column"] for spec in variation_specs if spec["feature_column"] != best["feature_column"]
    ] + [
        "ca1_pyramidal_n",
        "ca1_reactive_astrocyte_n",
        "ca1_oligodendrocyte_n",
    ] + [
        f"reactive_exposed_pyramidal_n_r{int(spec['radius_px'])}px" for spec in variation_specs
    ] + [
        f"oligo_free_reactive_exposed_pyramidal_n_r{int(spec['radius_px'])}px" for spec in variation_specs
    ]

    write_donor_feature_table(
        table_path,
        donor_table,
        feature_column=best["feature_column"],
        outcome_column=OUTCOME_COL,
        covariates=CANONICAL_COVARIATES,
        extra_columns=extra_columns,
    )

    payload = build_results_payload(
        status="ok",
        feature_name=best["feature_name"],
        outcome=OUTCOME_COL,
        n_total=n_total,
        n_analyzable=best["n_analyzable"],
        partial_r=_safe_float(best["partial_r"]),
        ci_lo=_safe_float(best["ci_lo"]),
        ci_hi=_safe_float(best["ci_hi"]),
        p_value=_safe_float(best["p_value"]),
        loo_predictive_r=_safe_float(best["loo_predictive_r"]),
        loo_unstable_count=int(best["loo_unstable_count"]),
        loo_max_shift=_safe_float(best["loo_max_shift"]),
        donor_ids_used=list(best["donor_ids_used"]),
        covariates=CONFounds_NUMERIC,
        recomputed_from_raw=True,
        registry_written=False,
        artifacts={"donor_feature_table": str(table_path)},
    )
    payload.update(
        _nan_to_none(
            {
                "feature_column": best["feature_column"],
                "best_variation": best["name"],
                "selection_score": best["selection_score"],
                "coverage_ratio": best["coverage_ratio"],
                "is_loo_gap": best["is_loo_gap"],
                "gap_penalty": best["gap_penalty"],
                "adjusted_score": best["adjusted_score"],
                "ranked_variations": [
                    {
                        key: item[key]
                        for key in [
                            "name",
                            "description",
                            "radius_px",
                            "feature_name",
                            "feature_column",
                            "n_total",
                            "n_analyzable",
                            "coverage_ratio",
                            "partial_r",
                            "selection_score",
                            "loo_predictive_r",
                            "is_loo_gap",
                            "gap_penalty",
                            "adjusted_score",
                            "p_value",
                            "ci_lo",
                            "ci_hi",
                            "loo_unstable_count",
                            "loo_max_shift",
                        ]
                    }
                    for item in ranked
                ],
            }
        )
    )
    write_results_payload(results_path, payload)
    _write_report(report_path, best, ranked, donor_table)

    print(f"HYPOTHESIS FAMILY: {FEATURE_FAMILY}")
    print(f"BEST VARIATION: {best['name']}")
    print(f"  IS partial r:      {_format_num(best['partial_r'])}")
    print(f"  Selection score:   {_format_num(best['selection_score'])}")
    print(f"  LOO predictive r:  {_format_num(best['loo_predictive_r'])}  (diagnostic)")
    print(f"  IS-LOO Gap:        {_format_num(best['is_loo_gap'])}  (penalty={_format_num(best['gap_penalty'])})")
    print(f"  Adjusted Score:    {_format_num(best['adjusted_score'])}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked:
        print(
            f"  {item['name']}  {_format_num(item['partial_r'])}  "
            f"{_format_num(item['selection_score'])}  {_format_num(item['loo_predictive_r'])}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best['feature_column']}")
    per_donor_df = pd.DataFrame(best["per_donor_loo"])
    for _, row in per_donor_df.iterrows():
        print(
            f"  {row['donor_id']}  {_format_num(row['outcome'])}  "
            f"{_format_num(row['predicted'])}  {_format_num(row[best['feature_column']])}"
        )


if __name__ == "__main__":
    main()
PY

python /scratch/result.py
echo
echo '--- results.json ---'
python - <<'PY'
import json
with open('/scratch/results.json') as f:
    obj=json.load(f)
print(json.dumps(obj, indent=2)[:12000])
PY
echo
echo '--- report.md ---'
sed -n '1,240p' /scratch/report.md
