set -e
cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import sys
import warnings
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, "/shared/lib")

from shared_analysis import (  # type: ignore
    bootstrap_partial_correlation,
    build_cell_table,
    leave_one_out_summary,
    load_training_cohort,
    write_donor_feature_table,
)
from shared_analysis.stats import (  # type: ignore
    bootstrap_partial_correlation_stability,
    partial_correlation,
    residualized_loo_predictive_correlation,
)

warnings.filterwarnings("ignore", message="Object at .* is not recognized as a component of a Zarr hierarchy\\.")
warnings.filterwarnings("ignore", category=UserWarning)

HERE = Path(__file__).resolve().parent
DATA_ROOT = Path("/data")
RESULTS_PATH = HERE / "results.json"
REPORT_PATH = HERE / "report.md"
DONOR_TABLE_PATH = HERE / "donor_feature_table.csv"

OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]

FAMILY_NAME = "ca1_reactive_hotspot_lymphocyte_enrichment"
BASELINE_VARIATION = "reactive_minus_astro_hotspot_local_lymphocyte_fraction_r70px_k3"
CANONICAL_FALLBACK_VARIATION = BASELINE_VARIATION

HOTSPOT_RADIUS_PX = 70.0
HOTSPOT_MIN_NEIGHBORS = 3

VARIATIONS: dict[str, dict[str, Any]] = {
    "reactive_minus_astro_hotspot_local_lymphocyte_fraction_r70px_k3": {
        "radius_px": 70.0,
        "k_same_class": 3,
        "feature_column": (
            "ca1_reactive_hotspot_lymphocyte_enrichment__"
            "reactive_minus_astro_hotspot_local_lymphocyte_fraction_r70px_k3"
        ),
        "description": (
            "Mean local lymphocyte fraction within 70 px around CA1 reactive-astrocyte hotspots "
            "minus the same fraction around matched CA1 astrocyte hotspots; hotspot rule uses "
            ">=3 same-class neighbors within 70 px."
        ),
    },
    "reactive_minus_astro_hotspot_local_lymphocyte_fraction_r90px_k3": {
        "radius_px": 90.0,
        "k_same_class": 3,
        "feature_column": (
            "ca1_reactive_hotspot_lymphocyte_enrichment__"
            "reactive_minus_astro_hotspot_local_lymphocyte_fraction_r90px_k3"
        ),
        "description": (
            "Mean local lymphocyte fraction within 90 px around CA1 reactive-astrocyte hotspots "
            "minus the same fraction around matched CA1 astrocyte hotspots; hotspot rule uses "
            ">=3 same-class neighbors within 70 px."
        ),
    },
}


@lru_cache(maxsize=1)
def _training_cohort(data_root: str | Path = DATA_ROOT) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = (
        cohort["sex"].astype(str).str.strip().str.lower().map({"female": 0.0, "male": 1.0}).astype(float)
    )
    return cohort


def _safe_float(value: Any) -> float:
    try:
        value = float(value)
    except Exception:
        return float("nan")
    return value if math.isfinite(value) else float("nan")


def _fmt(value: Any, digits: int = 4) -> str:
    value = _safe_float(value)
    return f"{value:.{digits}f}" if math.isfinite(value) else "nan"


def _query_counts(tree: cKDTree | None, points: np.ndarray, radius: float) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    if tree is None or len(points) == 0:
        return np.zeros(len(points), dtype=float)
    try:
        counts = tree.query_ball_point(points, float(radius), return_length=True)
        return np.asarray(counts, dtype=float)
    except TypeError:
        neighborhoods = tree.query_ball_point(points, float(radius))
        return np.asarray([len(xs) for xs in neighborhoods], dtype=float)


def _mean_local_fraction(
    *,
    hotspot_points: np.ndarray,
    all_tree: cKDTree,
    target_tree: cKDTree | None,
    radius_px: float,
) -> float:
    hotspot_points = np.asarray(hotspot_points, dtype=float)
    if hotspot_points.ndim != 2 or len(hotspot_points) == 0:
        return float("nan")
    denom = _query_counts(all_tree, hotspot_points, radius_px)
    numer = _query_counts(target_tree, hotspot_points, radius_px) if target_tree is not None else np.zeros(len(hotspot_points), dtype=float)
    mask = denom > 0
    if not np.any(mask):
        return float("nan")
    fractions = numer[mask] / denom[mask]
    return float(np.mean(fractions)) if len(fractions) else float("nan")


def _compute_slide_variation_scores(slide_path: str | Path) -> dict[str, Any]:
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    needed = ["x", "y", "cell_type", "region"]
    cells = cells.loc[:, needed].copy()
    cells = cells.loc[(cells["region"] == "CA1") & cells["cell_type"].notna()].reset_index(drop=True)

    result: dict[str, Any] = {
        "n_ca1_cells": int(len(cells)),
        "n_ca1_lymphocyte": 0,
        "n_ca1_reactive": 0,
        "n_ca1_astro": 0,
        "n_reactive_hotspots": 0,
        "n_astro_hotspots": 0,
        "variation_scores": {name: float("nan") for name in VARIATIONS},
    }
    if cells.empty:
        return result

    points_all = cells[["x", "y"]].to_numpy(dtype=np.float32, copy=False)
    cell_types = cells["cell_type"].astype(str).to_numpy()

    reactive_points = points_all[cell_types == "Reactive Astrocyte"]
    astro_points = points_all[cell_types == "Astrocyte"]
    lymph_points = points_all[cell_types == "Lymphocyte"]

    result["n_ca1_lymphocyte"] = int(len(lymph_points))
    result["n_ca1_reactive"] = int(len(reactive_points))
    result["n_ca1_astro"] = int(len(astro_points))

    if len(reactive_points) == 0 or len(astro_points) == 0:
        return result

    all_tree = cKDTree(points_all)
    lymph_tree = cKDTree(lymph_points) if len(lymph_points) else None

    reactive_tree = cKDTree(reactive_points)
    astro_tree = cKDTree(astro_points)

    reactive_same_counts = _query_counts(reactive_tree, reactive_points, HOTSPOT_RADIUS_PX) - 1.0
    astro_same_counts = _query_counts(astro_tree, astro_points, HOTSPOT_RADIUS_PX) - 1.0

    reactive_hotspots = reactive_points[reactive_same_counts >= HOTSPOT_MIN_NEIGHBORS]
    astro_hotspots = astro_points[astro_same_counts >= HOTSPOT_MIN_NEIGHBORS]

    result["n_reactive_hotspots"] = int(len(reactive_hotspots))
    result["n_astro_hotspots"] = int(len(astro_hotspots))

    if len(reactive_hotspots) == 0 or len(astro_hotspots) == 0:
        return result

    scores: dict[str, float] = {}
    for variation_name, spec in VARIATIONS.items():
        radius_px = float(spec["radius_px"])
        reactive_mean = _mean_local_fraction(
            hotspot_points=reactive_hotspots,
            all_tree=all_tree,
            target_tree=lymph_tree,
            radius_px=radius_px,
        )
        astro_mean = _mean_local_fraction(
            hotspot_points=astro_hotspots,
            all_tree=all_tree,
            target_tree=lymph_tree,
            radius_px=radius_px,
        )
        if not math.isfinite(reactive_mean) or not math.isfinite(astro_mean):
            scores[variation_name] = float("nan")
        else:
            scores[variation_name] = float(reactive_mean - astro_mean)
    result["variation_scores"] = scores
    return result


def _compute_all_donor_features(data_root: Path) -> pd.DataFrame:
    cohort = _training_cohort(data_root).copy()
    rows: list[dict[str, Any]] = []
    for row in cohort.itertuples(index=False):
        donor_id = str(row.donor_id)
        slide_name = str(row.slide_name)
        slide_result = _compute_slide_variation_scores(data_root / slide_name)

        out: dict[str, Any] = {
            "donor_id": donor_id,
            "slide_name": slide_name,
            OUTCOME_COLUMN: _safe_float(getattr(row, OUTCOME_COLUMN)),
            "max_age_vis": _safe_float(getattr(row, "max_age_vis")),
            "braak_numeric": _safe_float(getattr(row, "braak_numeric")),
            "cerad_ordinal": _safe_float(getattr(row, "cerad_ordinal")),
            "sex_binary": _safe_float(getattr(row, "sex_binary")),
            "sex": str(getattr(row, "sex")),
            "cognitive_status": str(getattr(row, "cognitive_status")),
            "n_ca1_cells": int(slide_result["n_ca1_cells"]),
            "n_ca1_lymphocyte": int(slide_result["n_ca1_lymphocyte"]),
            "n_ca1_reactive": int(slide_result["n_ca1_reactive"]),
            "n_ca1_astro": int(slide_result["n_ca1_astro"]),
            "n_reactive_hotspots": int(slide_result["n_reactive_hotspots"]),
            "n_astro_hotspots": int(slide_result["n_astro_hotspots"]),
        }
        for variation_name, spec in VARIATIONS.items():
            out[spec["feature_column"]] = _safe_float(slide_result["variation_scores"].get(variation_name))
        rows.append(out)
    return pd.DataFrame(rows)


def _selection_score(partial_r: Any, n_analyzable: Any, n_total: Any) -> float:
    partial_r = _safe_float(partial_r)
    n_analyzable = _safe_float(n_analyzable)
    n_total = _safe_float(n_total)
    if not math.isfinite(partial_r) or not math.isfinite(n_analyzable) or not math.isfinite(n_total) or n_total <= 0:
        return float("nan")
    return float(abs(partial_r) * max(0.0, min(1.0, n_analyzable / n_total)))


def _adjusted_metrics(partial_r: Any, loo_predictive_r: Any) -> tuple[float, float, float]:
    partial_r = _safe_float(partial_r)
    loo_predictive_r = _safe_float(loo_predictive_r)
    if not math.isfinite(partial_r) or not math.isfinite(loo_predictive_r):
        return float("nan"), float("nan"), float("nan")
    gap = float(abs(partial_r) - abs(loo_predictive_r))
    penalty = float(max(0.0, gap - 0.15) * 0.5)
    adjusted = float(abs(loo_predictive_r) - penalty)
    return gap, penalty, adjusted


def _loo_prediction_table(frame: pd.DataFrame, *, feature_col: str) -> pd.DataFrame:
    cols = ["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]
    clean = frame.loc[:, cols].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    rows: list[dict[str, Any]] = []
    for idx in range(len(clean)):
        train = clean.drop(index=idx).reset_index(drop=True)
        test = clean.iloc[[idx]].reset_index(drop=True)

        x_train = np.column_stack([np.ones(len(train), dtype=float), train[CONFOUND_COLUMNS].to_numpy(dtype=float)])
        x_test = np.column_stack([np.ones(len(test), dtype=float), test[CONFOUND_COLUMNS].to_numpy(dtype=float)])

        beta_feature, *_ = np.linalg.lstsq(x_train, train[feature_col].to_numpy(dtype=float), rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train, train[OUTCOME_COLUMN].to_numpy(dtype=float), rcond=None)

        resid_feature_train = train[feature_col].to_numpy(dtype=float) - x_train @ beta_feature
        resid_outcome_train = train[OUTCOME_COLUMN].to_numpy(dtype=float) - x_train @ beta_outcome
        denom = float(np.dot(resid_feature_train, resid_feature_train))
        slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom) if denom > 0 else 0.0

        conf_pred_outcome = float((x_test @ beta_outcome)[0])
        conf_pred_feature = float((x_test @ beta_feature)[0])
        feature_value = float(test.loc[0, feature_col])
        predicted = conf_pred_outcome + slope * (feature_value - conf_pred_feature)
        actual = float(test.loc[0, OUTCOME_COLUMN])

        rows.append(
            {
                "donor_id": str(test.loc[0, "donor_id"]),
                "outcome": actual,
                "predicted": float(predicted),
                feature_col: feature_value,
                "abs_error": abs(actual - float(predicted)),
            }
        )
    return pd.DataFrame(rows)


def _evaluate_variation(donor_table: pd.DataFrame, variation_name: str) -> dict[str, Any]:
    spec = VARIATIONS[variation_name]
    feature_col = str(spec["feature_column"])

    stats = partial_correlation(
        donor_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    boot = bootstrap_partial_correlation(
        donor_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        n_boot=1000,
        random_state=variation_name.__hash__() % (2**32),
    )
    stability = bootstrap_partial_correlation_stability(
        donor_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        n_boot=400,
        sample_frac=0.8,
        random_state=(len(variation_name) * 17) % (2**32),
    )
    loo_r = residualized_loo_predictive_correlation(
        donor_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    loo_summary = leave_one_out_summary(
        donor_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        id_col="donor_id",
        unstable_delta=0.10,
    )
    per_donor_loo = _loo_prediction_table(donor_table, feature_col=feature_col)

    n_total = int(len(donor_table))
    analyzable = donor_table[[feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]].replace([np.inf, -np.inf], np.nan).dropna()
    n_analyzable = int(len(analyzable))
    coverage_ratio = float(n_analyzable / n_total) if n_total > 0 else float("nan")
    selection = _selection_score(stats.get("partial_r"), n_analyzable, n_total)
    gap, penalty, adjusted = _adjusted_metrics(stats.get("partial_r"), loo_r)

    return {
        "name": variation_name,
        "feature_name": feature_col,
        "feature_column": feature_col,
        "description": spec["description"],
        "radius_px": float(spec["radius_px"]),
        "k_same_class": int(spec["k_same_class"]),
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "coverage_ratio": coverage_ratio,
        "partial_r": _safe_float(stats.get("partial_r")),
        "p_value": _safe_float(stats.get("p_value")),
        "ci_lo": _safe_float(boot.get("ci_lo")),
        "ci_hi": _safe_float(boot.get("ci_hi")),
        "selection_score": selection,
        "loo_predictive_r": _safe_float(loo_r),
        "is_loo_gap": gap,
        "gap_penalty": penalty,
        "adjusted_score": adjusted,
        "bootstrap_median_partial_r": _safe_float(stability.get("bootstrap_median_partial_r")),
        "bootstrap_sign_consistency": _safe_float(stability.get("bootstrap_sign_consistency")),
        "bootstrap_q25_partial_r": _safe_float(stability.get("bootstrap_q25_partial_r")),
        "bootstrap_q75_partial_r": _safe_float(stability.get("bootstrap_q75_partial_r")),
        "loo_unstable_count": int(loo_summary.get("unstable_count", 0) or 0),
        "loo_max_shift": _safe_float(loo_summary.get("max_shift")),
        "loo_unstable_donors": loo_summary.get("unstable_donors", []),
        "per_donor_loo": per_donor_loo.to_dict(orient="records"),
    }


def _ranking_key(item: dict[str, Any]) -> tuple[float, float, float]:
    sel = _safe_float(item.get("selection_score"))
    pr = abs(_safe_float(item.get("partial_r")))
    loo = abs(_safe_float(item.get("loo_predictive_r")))
    if not math.isfinite(sel):
        sel = -1.0
    if not math.isfinite(pr):
        pr = -1.0
    if not math.isfinite(loo):
        loo = -1.0
    return (sel, pr, loo)


def _canonical_variation_name() -> str:
    if RESULTS_PATH.exists():
        try:
            payload = json.loads(RESULTS_PATH.read_text(encoding="utf-8"))
            best = str(payload.get("best_variation") or "").strip()
            if best in VARIATIONS:
                return best
        except Exception:
            pass
    return CANONICAL_FALLBACK_VARIATION


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Replay the winning CA1 reactive-hotspot lymphocyte-enrichment biomarker.
    """
    data_root = Path(data_root)
    cohort = _training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"].astype(str) == str(donor_id)]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    best_variation = _canonical_variation_name()
    scores = _compute_slide_variation_scores(data_root / slide_name)["variation_scores"]
    value = scores.get(best_variation)
    return float(value) if value is not None else float("nan")


def _json_default(obj: Any):
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj)!r} is not JSON serializable")


def _write_report(best: dict[str, Any], ranked: list[dict[str, Any]], donor_table: pd.DataFrame) -> None:
    best_feature_col = str(best["feature_column"])
    loo = pd.DataFrame(best.get("per_donor_loo", []))
    worst = loo.sort_values("abs_error", ascending=False).head(3) if not loo.empty else pd.DataFrame()

    other_lines = []
    for item in ranked[1:]:
        other_lines.append(
            f"- `{item['name']}`: partial_r={_fmt(item['partial_r'])}, selection_score={_fmt(item['selection_score'])}, "
            f"loo_predictive_r={_fmt(item['loo_predictive_r'])}"
        )
    if not other_lines:
        other_lines.append("- No other planned variations were evaluated.")

    if worst.empty:
        error_text = "No analyzable LOO error table was available."
        next_text = "Repeat the same niche family only after confirming donor analyzability."
    else:
        donors = ", ".join(str(x) for x in worst["donor_id"].tolist())
        merged = donor_table.merge(loo[["donor_id", "predicted", "abs_error"]], on="donor_id", how="left")
        top = merged.set_index("donor_id").loc[worst["donor_id"].tolist()].reset_index()
        overcalled = top.loc[top["predicted"] < top[OUTCOME_COLUMN]]
        undercalled = top.loc[top["predicted"] > top[OUTCOME_COLUMN]]
        shared_bits: list[str] = []
        if len(overcalled) > 0:
            shared_bits.append("several were more impaired than predicted, suggesting decline not fully captured by focal lymphocyte enrichment")
        if len(undercalled) > 0:
            shared_bits.append("others carried the niche signal without equally severe decline, suggesting inflammatory-context signal can decouple from clinical slope")
        shared = "; ".join(shared_bits) if shared_bits else "the main shared feature was simply high absolute residual error"
        error_text = f"The largest LOO errors were {donors}; {shared}."
        next_text = (
            "Keep the same CA1 reactive-hotspot framework, but next sweep whether the immune context is better captured by a tighter "
            "cell-state gate inside the same family, for example hotspot-local lymphocyte presence versus fraction or a reactive-hotspot "
            "subset with stronger same-class crowding."
        )

    winner_name = str(best["name"])
    loser_name = str(ranked[1]["name"]) if len(ranked) > 1 else "the nearby alternative"
    if winner_name.endswith("r70px_k3"):
        why_worked = "The tighter 70 px readout preserved the most focal immune niche contrast around reactive hotspots."
        why_failed = f"The broader {loser_name} window likely diluted a sparse lymphocyte signal by averaging in surrounding CA1 background."
    else:
        why_worked = "The broader 90 px readout stabilized a sparse lymphocyte context better than the tighter neighborhood."
        why_failed = f"The tighter {loser_name} window was likely too noisy for a rare-cell immune readout."

    report = f"""## Summary
One sentence: tested the `{FAMILY_NAME}` family, the winning variation was `{winner_name}`, and the local winner selection score was `{_fmt(best['selection_score'])}`.

## Metrics
Winning variation `{winner_name}` (`{best_feature_col}`): IS partial r = {_fmt(best['partial_r'])}, selection score = {_fmt(best['selection_score'])}, LOO predictive r = {_fmt(best['loo_predictive_r'])}, IS-LOO gap = {_fmt(best['is_loo_gap'])} (penalty={_fmt(best['gap_penalty'])}), adjusted score = {_fmt(best['adjusted_score'])}.
Other ranked variations:
{chr(10).join(other_lines)}

## Findings
1. What worked and why (tie to the biological meaning of the target): {why_worked} Biologically, the feature asks whether CA1 reactive astrocyte microclusters sit in a more lymphocyte-rich microenvironment than matched homeostatic astrocyte microclusters.
2. What failed and why (specific to the chosen hypothesis and what went wrong): {why_failed} The losing local sweep stayed in the same niche family but changed only the measurement radius.
3. Error pattern: {error_text}

## Rationale
The best variation is biologically coherent because it isolates a specific inflammatory niche: hotspot-positive Reactive Astrocytes in CA1, benchmarked against hotspot-positive Astrocytes from the same region, and summarized as a donor-level delta in local lymphocyte fraction. That comparison controls away global cell-density and region-size effects better than a whole-CA1 lymphocyte summary. It beat the nearby alternative because its radius gave the better balance between focal specificity and stability for a rare immune-cell context. It may add information beyond the current panel because prior kept features in this hotspot series emphasized neuronal loss, oligodendrocyte depletion, and corpora-amylacea context rather than explicit lymphocyte enrichment.

## Interpretation
The signal seems to mean that some donors have CA1 reactive astrocyte hotspots embedded in a more immune-infiltrated local niche than ordinary astrocyte hotspots. Population: Reactive Astrocyte versus Astrocyte hotspot-positive cells in CA1. Niche: local CA1 neighborhoods within the winning radius around those hotspot centers. Feature summary: donor-level mean local lymphocyte fraction around reactive hotspots minus the same mean around astrocyte hotspots. Simplest observable pattern: focal reactive astrocyte microclusters that appear to sit nearer small lymphocyte pockets than nearby baseline astrocyte clusters.

## Next
{next_text}
"""
    REPORT_PATH.write_text(report, encoding="utf-8")


def main() -> int:
    donor_table = _compute_all_donor_features(DATA_ROOT)

    ranked = [_evaluate_variation(donor_table, name) for name in VARIATIONS]
    ranked.sort(key=_ranking_key, reverse=True)
    best = ranked[0]
    best_feature_col = str(best["feature_column"])
    best_variation = str(best["name"])

    donor_feature_table = donor_table.copy()
    write_donor_feature_table(
        DONOR_TABLE_PATH,
        donor_feature_table,
        feature_column=best_feature_col,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUND_COLUMNS,
        extra_columns=[
            "n_ca1_cells",
            "n_ca1_lymphocyte",
            "n_ca1_reactive",
            "n_ca1_astro",
            "n_reactive_hotspots",
            "n_astro_hotspots",
        ],
    )

    payload = {
        "status": "ok",
        "feature_name": best_feature_col,
        "family_name": FAMILY_NAME,
        "outcome": OUTCOME_COLUMN,
        "covariates": CONFOUND_COLUMNS,
        "best_variation": best_variation,
        "feature_column": best_feature_col,
        "n_total": int(best["n_total"]),
        "n_analyzable": int(best["n_analyzable"]),
        "coverage_ratio": _safe_float(best["coverage_ratio"]),
        "partial_r": _safe_float(best["partial_r"]),
        "p_value": _safe_float(best["p_value"]),
        "ci_lo": _safe_float(best["ci_lo"]),
        "ci_hi": _safe_float(best["ci_hi"]),
        "selection_score": _safe_float(best["selection_score"]),
        "loo_predictive_r": _safe_float(best["loo_predictive_r"]),
        "is_loo_gap": _safe_float(best["is_loo_gap"]),
        "gap_penalty": _safe_float(best["gap_penalty"]),
        "adjusted_score": _safe_float(best["adjusted_score"]),
        "bootstrap_median_partial_r": _safe_float(best["bootstrap_median_partial_r"]),
        "bootstrap_sign_consistency": _safe_float(best["bootstrap_sign_consistency"]),
        "loo_unstable_count": int(best["loo_unstable_count"]),
        "loo_max_shift": _safe_float(best["loo_max_shift"]),
        "ranked_variations": ranked,
        "artifacts": {
            "donor_feature_table": str(DONOR_TABLE_PATH),
            "report": str(REPORT_PATH),
            "result_py": str(HERE / "result.py"),
        },
    }
    RESULTS_PATH.write_text(json.dumps(payload, indent=2, default=_json_default) + "\n", encoding="utf-8")

    _write_report(best, ranked, donor_table)

    print(f"HYPOTHESIS FAMILY: {FAMILY_NAME}")
    print(f"BEST VARIATION: {best_variation}")
    print(f"  IS partial r:      {_fmt(best['partial_r'])}")
    print(f"  Selection score:   {_fmt(best['selection_score'])}")
    print(f"  LOO predictive r:  {_fmt(best['loo_predictive_r'])}  (diagnostic)")
    print(f"  IS-LOO Gap:        {_fmt(best['is_loo_gap'])}  (penalty={_fmt(best['gap_penalty'])})")
    print(f"  Adjusted Score:    {_fmt(best['adjusted_score'])}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked:
        print(
            f"  {item['name']}  {_fmt(item['partial_r'])}  {_fmt(item['selection_score'])}  {_fmt(item['loo_predictive_r'])}"
        )
    print()
    print("PER-DONOR (LOO):")
    loo = pd.DataFrame(best.get("per_donor_loo", []))
    if not loo.empty:
        for row in loo.itertuples(index=False):
            rowd = row._asdict()
            print(
                f"  {rowd['donor_id']}  {_fmt(rowd['outcome'])}  {_fmt(rowd['predicted'])}  {_fmt(rowd[best_feature_col])}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
PY
python /scratch/result.py
echo
echo '--- results.json ---'
python - <<'PY'
import json
obj=json.load(open('/scratch/results.json'))
print(json.dumps({k: obj[k] for k in ['best_variation','feature_column','selection_score','partial_r','loo_predictive_r','n_analyzable','n_total']}, indent=2))
print('ranked_variations:', [x['name'] for x in obj['ranked_variations']])
PY
echo
echo '--- report.md head ---'
sed -n '1,220p' /scratch/report.md
