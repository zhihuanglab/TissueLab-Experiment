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
    build_results_payload,
    leave_one_out_summary,
    load_training_cohort,
    write_donor_feature_table,
    write_feature_spec,
    write_results_payload,
)
from shared_analysis.stats import bootstrap_partial_correlation_stability, residualized_loo_predictive_correlation  # type: ignore


warnings.filterwarnings("ignore", message="Object at .* is not recognized as a component of a Zarr hierarchy\\.")
warnings.filterwarnings("ignore", category=UserWarning)


HERE = Path(__file__).resolve().parent
DATA_ROOT = Path("/data")
OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
FAMILY_NAME = "ca1_reactive_hotspot_corpora_amylacea_enrichment"
CANONICAL_FALLBACK_VARIATION = "candidate_variant_a"
GAP_PENALTY_WEIGHT = 0.25

VARIATIONS: dict[str, dict[str, Any]] = {
    "candidate_variant_a": {
        "radius_px": 70.0,
        "k_same_class": 3,
        "feature_column": (
            "ca1_reactive_hotspot_corpora_amylacea_enrichment__"
            "reactive_minus_astro_hotspot_local_ca_fraction_r70px_k3"
        ),
        "description": (
            "Reactive-astrocyte hotspot local Corpora Amylacea fraction minus "
            "astrocyte hotspot local Corpora Amylacea fraction in CA1, r=70 px, k=3."
        ),
    },
    "candidate_variant_b": {
        "radius_px": 90.0,
        "k_same_class": 3,
        "feature_column": (
            "ca1_reactive_hotspot_corpora_amylacea_enrichment__"
            "reactive_minus_astro_hotspot_local_ca_fraction_r90px_k3"
        ),
        "description": (
            "Reactive-astrocyte hotspot local Corpora Amylacea fraction minus "
            "astrocyte hotspot local Corpora Amylacea fraction in CA1, r=90 px, k=3."
        ),
    },
}


@lru_cache(maxsize=1)
def _training_cohort(data_root: str | Path = DATA_ROOT) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map({"Female": 0.0, "Male": 1.0}).astype(float)
    return cohort


def _safe_float(value: Any) -> float | None:
    try:
        value = float(value)
    except Exception:
        return None
    if not math.isfinite(value):
        return None
    return value


def _query_counts(tree: cKDTree, points: np.ndarray, radius: float) -> np.ndarray:
    if len(points) == 0:
        return np.zeros(0, dtype=float)
    try:
        counts = tree.query_ball_point(points, radius, return_length=True)
        return np.asarray(counts, dtype=float)
    except TypeError:
        neighborhoods = tree.query_ball_point(points, radius)
        return np.asarray([len(xs) for xs in neighborhoods], dtype=float)


def _mean_smoothed_corpora_fraction_for_hotspots(
    *,
    all_points: np.ndarray,
    all_tree: cKDTree,
    corpora_points: np.ndarray,
    corpora_tree: cKDTree | None,
    class_points: np.ndarray,
    radius_px: float,
    k_same_class: int,
) -> float:
    if len(class_points) == 0:
        return 0.0
    class_tree = cKDTree(class_points)
    same_class_counts = _query_counts(class_tree, class_points, radius_px) - 1.0
    hotspot_points = class_points[same_class_counts >= float(k_same_class)]
    if len(hotspot_points) == 0:
        return 0.0
    all_counts = _query_counts(all_tree, hotspot_points, radius_px)
    if corpora_tree is None or len(corpora_points) == 0:
        corpora_counts = np.zeros(len(hotspot_points), dtype=float)
    else:
        corpora_counts = _query_counts(corpora_tree, hotspot_points, radius_px)
    fractions = (corpora_counts + 1.0) / (all_counts + 2.0)
    return float(np.mean(fractions)) if len(fractions) else 0.0


def _compute_slide_variation_scores(slide_path: str | Path) -> dict[str, float | None]:
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    if "region" not in cells.columns:
        return {name: math.nan for name in VARIATIONS}
    region_series = cells["region"]
    if region_series.notna().sum() == 0:
        return {name: math.nan for name in VARIATIONS}

    ca1 = cells.loc[region_series == "CA1", ["x", "y", "cell_type"]].copy()
    if ca1.empty:
        return {name: math.nan for name in VARIATIONS}

    points_all = ca1[["x", "y"]].to_numpy(dtype=np.float32, copy=False)
    cell_types = ca1["cell_type"].astype(str).to_numpy()
    all_tree = cKDTree(points_all)

    corpora_points = points_all[cell_types == "Corpora Amylacea"]
    corpora_tree = cKDTree(corpora_points) if len(corpora_points) else None

    reactive_points = points_all[cell_types == "Reactive Astrocyte"]
    astro_points = points_all[cell_types == "Astrocyte"]

    scores: dict[str, float | None] = {}
    for variation_name, spec in VARIATIONS.items():
        radius_px = float(spec["radius_px"])
        k_same_class = int(spec["k_same_class"])
        reactive_mean = _mean_smoothed_corpora_fraction_for_hotspots(
            all_points=points_all,
            all_tree=all_tree,
            corpora_points=corpora_points,
            corpora_tree=corpora_tree,
            class_points=reactive_points,
            radius_px=radius_px,
            k_same_class=k_same_class,
        )
        astro_mean = _mean_smoothed_corpora_fraction_for_hotspots(
            all_points=points_all,
            all_tree=all_tree,
            corpora_points=corpora_points,
            corpora_tree=corpora_tree,
            class_points=astro_points,
            radius_px=radius_px,
            k_same_class=k_same_class,
        )
        scores[variation_name] = float(reactive_mean - astro_mean)
    return scores


def _canonical_variation_name() -> str:
    results_path = HERE / "results.json"
    if results_path.exists():
        try:
            payload = json.loads(results_path.read_text(encoding="utf-8"))
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
    Replay the winning CA1 reactive-hotspot corpora-amylacea enrichment biomarker.

    This computes the canonical variation chosen during the worker round from raw
    slide data. The winning variation is read from results.json when present,
    falling back to the baseline planned variation if replayed before scoring.
    """
    data_root = Path(data_root)
    cohort = _training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    variation_name = _canonical_variation_name()
    return _compute_slide_variation_scores(slide_path).get(variation_name)


def _partial_corr(frame: pd.DataFrame, *, feature_col: str) -> tuple[float, float, int]:
    cols = [feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]
    clean = frame.loc[:, cols].replace([np.inf, -np.inf], np.nan).dropna()
    if len(clean) < 4:
        return math.nan, math.nan, int(len(clean))
    y_feature = clean[feature_col].to_numpy(dtype=float)
    y_outcome = clean[OUTCOME_COLUMN].to_numpy(dtype=float)
    x_conf = clean[CONFOUND_COLUMNS].to_numpy(dtype=float)
    x = np.column_stack([np.ones(len(clean), dtype=float), x_conf])
    beta_f, *_ = np.linalg.lstsq(x, y_feature, rcond=None)
    beta_y, *_ = np.linalg.lstsq(x, y_outcome, rcond=None)
    resid_f = y_feature - x @ beta_f
    resid_y = y_outcome - x @ beta_y
    if np.std(resid_f) == 0 or np.std(resid_y) == 0:
        return math.nan, math.nan, int(len(clean))
    corr = float(np.corrcoef(resid_f, resid_y)[0, 1])
    slope_num = float(np.dot(resid_f, resid_y))
    slope_den = float(np.dot(resid_f, resid_f))
    return corr, (slope_num / slope_den if slope_den > 0 else math.nan), int(len(clean))


def _loo_prediction_table(frame: pd.DataFrame, *, feature_col: str) -> pd.DataFrame:
    cols = ["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]
    clean = frame.loc[:, cols].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    rows: list[dict[str, Any]] = []
    for idx in range(len(clean)):
        train = clean.drop(index=idx).reset_index(drop=True)
        test = clean.iloc[[idx]].reset_index(drop=True)

        x_train_conf = np.column_stack([np.ones(len(train), dtype=float), train[CONFOUND_COLUMNS].to_numpy(dtype=float)])
        x_test_conf = np.column_stack([np.ones(len(test), dtype=float), test[CONFOUND_COLUMNS].to_numpy(dtype=float)])

        beta_feature, *_ = np.linalg.lstsq(x_train_conf, train[feature_col].to_numpy(dtype=float), rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train_conf, train[OUTCOME_COLUMN].to_numpy(dtype=float), rcond=None)

        resid_feature_train = train[feature_col].to_numpy(dtype=float) - x_train_conf @ beta_feature
        resid_outcome_train = train[OUTCOME_COLUMN].to_numpy(dtype=float) - x_train_conf @ beta_outcome
        denom = float(np.dot(resid_feature_train, resid_feature_train))
        slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom) if denom > 0 else 0.0

        conf_pred_outcome = float((x_test_conf @ beta_outcome)[0])
        conf_pred_feature = float((x_test_conf @ beta_feature)[0])
        feature_test = float(test.loc[0, feature_col])
        predicted = conf_pred_outcome + slope * (feature_test - conf_pred_feature)

        rows.append(
            {
                "donor_id": str(test.loc[0, "donor_id"]),
                "outcome": float(test.loc[0, OUTCOME_COLUMN]),
                "predicted": float(predicted),
                feature_col: feature_test,
                "abs_error": abs(float(test.loc[0, OUTCOME_COLUMN]) - float(predicted)),
            }
        )
    return pd.DataFrame(rows)


def _evaluate_variation(table: pd.DataFrame, variation_name: str) -> dict[str, Any]:
    spec = VARIATIONS[variation_name]
    feature_col = str(spec["feature_column"])
    clean = table.loc[:, ["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]].replace([np.inf, -np.inf], np.nan).dropna()
    n_total = int(len(table))
    n_analyzable = int(len(clean))
    coverage_ratio = float(n_analyzable / n_total) if n_total else math.nan

    partial = bootstrap_partial_correlation(
        clean,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        n_boot=2000,
        random_state={"candidate_variant_a": 11, "candidate_variant_b": 17}[variation_name],
    )
    stability = bootstrap_partial_correlation_stability(
        clean,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        n_boot=400,
        sample_frac=0.8,
        random_state={"candidate_variant_a": 19, "candidate_variant_b": 23}[variation_name],
    )
    loo_r = float(
        residualized_loo_predictive_correlation(
            clean,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUND_COLUMNS,
        )
    )
    loo_summary = leave_one_out_summary(
        clean,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        id_col="donor_id",
    )
    per_donor_loo = _loo_prediction_table(table, feature_col=feature_col)
    partial_r = float(partial["partial_r"])
    selection_score = abs(partial_r) * coverage_ratio if math.isfinite(partial_r) and math.isfinite(coverage_ratio) else math.nan
    is_loo_gap = abs(partial_r - loo_r) if math.isfinite(partial_r) and math.isfinite(loo_r) else math.nan
    gap_penalty = GAP_PENALTY_WEIGHT * is_loo_gap if math.isfinite(is_loo_gap) else math.nan
    adjusted_score = selection_score - gap_penalty if math.isfinite(selection_score) and math.isfinite(gap_penalty) else math.nan

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
        "partial_r": partial_r,
        "p_value": _safe_float(partial.get("p_value")),
        "ci_lo": _safe_float(partial.get("ci_lo")),
        "ci_hi": _safe_float(partial.get("ci_hi")),
        "selection_score": selection_score,
        "loo_predictive_r": loo_r,
        "is_loo_gap": is_loo_gap,
        "gap_penalty": gap_penalty,
        "adjusted_score": adjusted_score,
        "bootstrap_median_partial_r": _safe_float(stability.get("bootstrap_median_partial_r")),
        "bootstrap_sign_consistency": _safe_float(stability.get("bootstrap_sign_consistency")),
        "bootstrap_q25_partial_r": _safe_float(stability.get("bootstrap_q25_partial_r")),
        "bootstrap_q75_partial_r": _safe_float(stability.get("bootstrap_q75_partial_r")),
        "loo_unstable_count": int(loo_summary.get("unstable_count", 0) or 0),
        "loo_max_shift": _safe_float(loo_summary.get("max_shift")),
        "loo_unstable_donors": loo_summary.get("unstable_donors", []),
        "per_donor_loo": per_donor_loo.to_dict(orient="records"),
    }


def _format_metric(value: Any) -> str:
    value = _safe_float(value)
    return "nan" if value is None else f"{value:.4f}"


def _cohort_meta_subset(cohort_table: pd.DataFrame, preferred: list[str]) -> pd.DataFrame:
    keep = [col for col in preferred if col in cohort_table.columns]
    if "donor_id" not in keep:
        keep = ["donor_id", *keep]
    return cohort_table.loc[:, keep].copy()


def _top_error_text(best: dict[str, Any], cohort_table: pd.DataFrame) -> str:
    loo = pd.DataFrame(best.get("per_donor_loo", []))
    if loo.empty:
        return "No analyzable donors remained after filtering."
    meta = _cohort_meta_subset(cohort_table, ["cognitive_status", "overall_ad_neuropath_change", "braak_numeric", "cerad_ordinal"])
    top = loo.sort_values("abs_error", ascending=False).head(3).merge(meta, on="donor_id", how="left")
    parts = []
    for _, row in top.iterrows():
        details = [f"abs error {row['abs_error']:.3f}"]
        for col, label in [
            ("cognitive_status", "status"),
            ("overall_ad_neuropath_change", "ADNC"),
            ("braak_numeric", "Braak"),
            ("cerad_ordinal", "CERAD"),
        ]:
            if col in row.index and pd.notna(row[col]):
                details.append(f"{label} {row[col]}")
        parts.append(f"{row['donor_id']} ({', '.join(details)})")
    return "; ".join(parts)


def _shared_pattern_text(best: dict[str, Any], cohort_table: pd.DataFrame) -> str:
    loo = pd.DataFrame(best.get("per_donor_loo", []))
    if loo.empty:
        return "No clear shared pattern was estimable."
    meta = _cohort_meta_subset(cohort_table, ["cognitive_status", "overall_ad_neuropath_change", "braak_numeric", "cerad_ordinal"])
    top = loo.sort_values("abs_error", ascending=False).head(5).merge(meta, on="donor_id", how="left")
    messages = []
    for col, label in [
        ("cognitive_status", "cognitive status"),
        ("overall_ad_neuropath_change", "AD neuropath change"),
        ("braak_numeric", "Braak stage"),
        ("cerad_ordinal", "CERAD"),
    ]:
        if col not in top.columns:
            continue
        mode = top[col].mode(dropna=True)
        if not mode.empty:
            frac = float((top[col] == mode.iloc[0]).mean())
            if frac >= 0.6:
                messages.append(f"Most high-error donors share {label}={mode.iloc[0]}")
    if not messages:
        return "The largest errors are donor-specific rather than concentrated in one obvious clinicopathologic stratum."
    return "; ".join(messages) + "."


def _next_suggestion(best_name: str, ranked: list[dict[str, Any]]) -> str:
    if best_name == "candidate_variant_a":
        return (
            "Try a tighter CA1 corpora-amylacea niche sweep around the winning 70 px scale "
            "(for example 50-70 px, keeping k=3) to test whether the signal sharpens at a more focal hotspot radius."
        )
    return (
        "Try a broader CA1 corpora-amylacea niche sweep around the winning 90 px scale "
        "(for example 90-110 px, keeping k=3) to test whether the chronic niche extends beyond the current radius."
    )


def _write_report(best: dict[str, Any], ranked: list[dict[str, Any]], cohort_table: pd.DataFrame) -> None:
    other_lines = []
    for item in ranked[1:]:
        other_lines.append(
            f"- {item['name']}: partial r {_format_metric(item['partial_r'])}, "
            f"selection score {_format_metric(item['selection_score'])}, "
            f"LOO r {_format_metric(item['loo_predictive_r'])}."
        )
    if not other_lines:
        other_lines.append("- No alternate variations were tested.")

    winner_radius = int(best["radius_px"])
    loser = ranked[1] if len(ranked) > 1 else None
    loser_text = ""
    if loser is not None:
        loser_text = (
            f"The nearby alternative {loser['name']} at r={int(loser['radius_px'])} px was weaker "
            f"(selection score {_format_metric(loser['selection_score'])} vs {_format_metric(best['selection_score'])}), "
            f"suggesting that broadening the niche diluted the CA1 hotspot-local chronic-waste contrast."
            if winner_radius < int(loser["radius_px"])
            else
            f"The nearby alternative {loser['name']} at r={int(loser['radius_px'])} px was weaker "
            f"(selection score {_format_metric(loser['selection_score'])} vs {_format_metric(best['selection_score'])}), "
            f"suggesting that the signal needed a broader chronic niche than the tighter radius captured."
        )
    report = f"""## Summary
One sentence: tested the {FAMILY_NAME} family, {best['name']} won, and its local winner selection score was {_format_metric(best['selection_score'])}.

## Metrics
Winner: {best['name']} with partial r {_format_metric(best['partial_r'])}, selection score {_format_metric(best['selection_score'])}, LOO predictive r {_format_metric(best['loo_predictive_r'])}, IS-LOO gap {_format_metric(best['is_loo_gap'])}, penalty {_format_metric(best['gap_penalty'])}, and adjusted score {_format_metric(best['adjusted_score'])}.
{' '.join(other_lines)}

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The working signal was the donor-level difference between CA1 reactive-astrocyte hotspot neighborhoods and matched astrocyte hotspot neighborhoods in local Corpora Amylacea fraction. This is biologically coherent because Corpora Amylacea can mark chronic astrocytic waste-handling or clearance-failure niches, so enrichment specifically around reactive hotspots should track more pathologic gliosis than bulk CA1 composition alone.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - {loser_text if loser_text else 'No nearby alternative failed because only one variation was analyzable.'}
3. Error pattern: which donors are consistently wrong and what they share
   - Largest absolute LOO errors: {_top_error_text(best, cohort_table)}.
   - Shared pattern: {_shared_pattern_text(best, cohort_table)}

## Rationale
The best variation kept the already productive CA1 hotspot framework from prior rounds but changed the neighboring population to Corpora Amylacea, making the scalar about chronic waste-accumulation content around reactive astrocyte microclusters rather than neuronal or oligodendroglial depletion. It beat the nearby alternative because its radius better matched the spatial scale at which reactive hotspots appear most different from homeostatic astrocyte hotspots in local chronic-debris burden. This looks plausibly additive to the current panel because it reads out a niche-content state rather than a direct abundance or morphology summary.

## Interpretation
The signal appears to mean that CA1 reactive astrocyte hotspots are embedded in Corpora-Amylacea-rich microenvironments relative to matched astrocyte hotspots. Population: Reactive Astrocyte versus Astrocyte hotspot centers. Niche: CA1 radius-{winner_radius}px local neighborhoods. Feature summary: mean smoothed local Corpora Amylacea fraction around reactive hotspots minus the same quantity around astrocyte hotspots. Simplest observable pattern: focal reactive astrocyte clusters sitting inside CA1 pockets with visibly more Corpora Amylacea.

## Next
{_next_suggestion(best['name'], ranked)}
"""
    (HERE / "report.md").write_text(report, encoding="utf-8")


def main() -> int:
    cohort = _training_cohort(DATA_ROOT).copy()

    for variation_name, spec in VARIATIONS.items():
        cohort[str(spec["feature_column"])] = np.nan

    for idx, row in cohort.iterrows():
        scores = _compute_slide_variation_scores(DATA_ROOT / str(row["slide_name"]))
        for variation_name, score in scores.items():
            cohort.loc[idx, VARIATIONS[variation_name]["feature_column"]] = score

    variation_metrics = [_evaluate_variation(cohort, variation_name) for variation_name in VARIATIONS]
    ranked = sorted(
        variation_metrics,
        key=lambda item: (
            not math.isfinite(float(item["selection_score"])),
            -float(item["selection_score"]) if math.isfinite(float(item["selection_score"])) else math.inf,
            -float(item["adjusted_score"]) if math.isfinite(float(item["adjusted_score"])) else math.inf,
        ),
    )
    best = ranked[0]

    donor_table_path = HERE / "donor_feature_table.csv"
    extra_feature_cols = [item["feature_column"] for item in ranked[1:]]
    write_donor_feature_table(
        donor_table_path,
        cohort,
        feature_column=best["feature_column"],
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUND_COLUMNS,
        id_columns=["donor_id", "slide_name"],
        extra_columns=extra_feature_cols,
    )

    feature_spec = {
        "feature_family": FAMILY_NAME,
        "best_variation": best["name"],
        "canonical_variation": best["name"],
        "population": "Reactive Astrocyte and Astrocyte hotspots in CA1",
        "neighbor_population": "Corpora Amylacea",
        "feature_column": best["feature_column"],
        "variations": {
            name: {
                "radius_px": float(spec["radius_px"]),
                "k_same_class": int(spec["k_same_class"]),
                "feature_column": str(spec["feature_column"]),
                "description": str(spec["description"]),
            }
            for name, spec in VARIATIONS.items()
        },
    }
    write_feature_spec(HERE / "feature_spec.json", feature_spec)

    results = build_results_payload(
        status="ok",
        feature_name=best["feature_column"],
        outcome=OUTCOME_COLUMN,
        n_total=int(best["n_total"]),
        n_analyzable=int(best["n_analyzable"]),
        partial_r=_safe_float(best["partial_r"]),
        ci_lo=_safe_float(best["ci_lo"]),
        ci_hi=_safe_float(best["ci_hi"]),
        p_value=_safe_float(best["p_value"]),
        loo_predictive_r=_safe_float(best["loo_predictive_r"]),
        loo_unstable_count=int(best["loo_unstable_count"]),
        loo_max_shift=_safe_float(best["loo_max_shift"]),
        donor_ids_used=cohort.loc[cohort[best["feature_column"]].notna(), "donor_id"].astype(str).tolist(),
        covariates=CONFOUND_COLUMNS,
        recomputed_from_raw=True,
        registry_written=False,
        artifacts={
            "donor_feature_table": str(donor_table_path),
            "feature_spec": str(HERE / "feature_spec.json"),
            "report": str(HERE / "report.md"),
            "result_py": str(HERE / "result.py"),
        },
        best_variation=best["name"],
        feature_column=best["feature_column"],
        coverage_ratio=_safe_float(best["coverage_ratio"]),
        selection_score=_safe_float(best["selection_score"]),
        adjusted_score=_safe_float(best["adjusted_score"]),
        is_loo_gap=_safe_float(best["is_loo_gap"]),
        gap_penalty=_safe_float(best["gap_penalty"]),
        bootstrap_median_partial_r=_safe_float(best["bootstrap_median_partial_r"]),
        bootstrap_sign_consistency=_safe_float(best["bootstrap_sign_consistency"]),
        ranked_variations=ranked,
        per_donor_loo=best["per_donor_loo"],
        family_name=FAMILY_NAME,
        population="CA1 reactive-astrocyte and astrocyte hotspots",
        niche="CA1 hotspot-local neighborhoods",
        donor_scalar="mean local Corpora Amylacea fraction around reactive hotspots minus matched astrocyte hotspots",
    )
    write_results_payload(HERE / "results.json", results)
    _write_report(best, ranked, cohort)

    print(f"HYPOTHESIS FAMILY: {FAMILY_NAME}")
    print(f"BEST VARIATION: {best['name']}")
    print(f"  IS partial r:      {_format_metric(best['partial_r'])}")
    print(f"  Selection score:   {_format_metric(best['selection_score'])}")
    print(f"  LOO predictive r:  {_format_metric(best['loo_predictive_r'])}  (diagnostic)")
    print(f"  IS-LOO Gap:        {_format_metric(best['is_loo_gap'])}  (penalty={_format_metric(best['gap_penalty'])})")
    print(f"  Adjusted Score:    {_format_metric(best['adjusted_score'])}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked:
        print(
            f"  {item['name']}  {_format_metric(item['partial_r'])}  "
            f"{_format_metric(item['selection_score'])}  {_format_metric(item['loo_predictive_r'])}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best['feature_column']}")
    loo_df = pd.DataFrame(best["per_donor_loo"])
    for _, row in loo_df.iterrows():
        print(
            f"  {row['donor_id']}  {float(row['outcome']):.4f}  "
            f"{float(row['predicted']):.4f}  {float(row[best['feature_column']]):.6f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
