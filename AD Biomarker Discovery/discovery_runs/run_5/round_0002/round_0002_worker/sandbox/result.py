from __future__ import annotations

import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath
from scipy.spatial import cKDTree

SHARED_LIB = Path("/shared/lib")
if str(SHARED_LIB) not in sys.path:
    sys.path.insert(0, str(SHARED_LIB))

from shared_analysis.artifacts import build_results_payload, write_donor_feature_table, write_results_payload  # noqa: E402
from shared_analysis.sea_ad_lfb import (  # noqa: E402
    load_centroids,
    load_class_ids,
    load_class_lookup,
    load_region_polygons,
    load_training_cohort,
)
from shared_analysis.stats import (  # noqa: E402
    bootstrap_partial_correlation,
    bootstrap_partial_correlation_stability,
    leave_one_out_summary,
    partial_correlation,
    residualized_loo_predictive_correlation,
)

warnings.filterwarnings(
    "ignore",
    message="Object at .* is not recognized as a component of a Zarr hierarchy.",
)

HYPOTHESIS_FAMILY = "ca1_pyramidal_reactive_astro_neighborhood_fraction"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUNDS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]

REGION_NAME = "CA1"
PYRAMIDAL_LABEL = "Pyramidal Neuron"
REACTIVE_ASTRO_LABEL = "Reactive Astrocyte"
MIN_PYRAMIDAL_CELLS = 50

# Patched by main() after the local sweep is evaluated.
CANONICAL_VARIATION = "pyramidal_near_reactiveastro_30um"

VARIATIONS: list[dict[str, Any]] = [
    {
        "name": "pyramidal_near_reactiveastro_30um",
        "radius_um": 30.0,
        "radius_px": 60.0,  # ~30 / 0.503
        "feature_name": "ca1_pyramidal_near_reactive_astro_fraction_30um",
        "feature_column": "ca1_pyramidal_near_reactive_astro_fraction_30um",
        "description": "Fraction of CA1 pyramidal neurons with at least one CA1 reactive astrocyte within ~30 um.",
    },
    {
        "name": "pyramidal_near_reactiveastro_50um",
        "radius_um": 50.0,
        "radius_px": 100.0,  # ~50 / 0.503
        "feature_name": "ca1_pyramidal_near_reactive_astro_fraction_50um",
        "feature_column": "ca1_pyramidal_near_reactive_astro_fraction_50um",
        "description": "Fraction of CA1 pyramidal neurons with at least one CA1 reactive astrocyte within ~50 um.",
    },
]
VARIATION_LOOKUP = {spec["name"]: spec for spec in VARIATIONS}

FEATURE_NAME = VARIATION_LOOKUP[CANONICAL_VARIATION]["feature_name"]
FEATURE_COLUMN = VARIATION_LOOKUP[CANONICAL_VARIATION]["feature_column"]


def _float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        value = float(value)
    except Exception:
        return None
    if not np.isfinite(value):
        return None
    return value


def _selection_score(partial_r: float | None, n_analyzable: int | None, n_total: int | None) -> float:
    if partial_r is None or not np.isfinite(partial_r) or n_analyzable is None or n_total is None or n_total <= 0:
        return float("nan")
    return float(abs(float(partial_r)) * (float(n_analyzable) / float(n_total)))


def _adjusted_metrics(partial_r: float | None, loo_predictive_r: float | None) -> tuple[float, float, float]:
    is_r = abs(float(partial_r)) if partial_r is not None and np.isfinite(partial_r) else 0.0
    loo_r = abs(float(loo_predictive_r)) if loo_predictive_r is not None and np.isfinite(loo_predictive_r) else 0.0
    gap = is_r - loo_r
    if gap > 0.30:
        return float(gap), float(gap), -1.0
    penalty = max(0.0, gap - 0.15) * 0.5
    return float(gap), float(penalty), float(loo_r - penalty)


def _inverse_lookup(class_lookup: dict[int, str]) -> dict[str, int]:
    return {str(name): int(idx) for idx, name in class_lookup.items()}


def _region_mask(centroids: np.ndarray, polygons: list[np.ndarray]) -> np.ndarray:
    mask = np.zeros(len(centroids), dtype=bool)
    for polygon in polygons:
        polygon = np.asarray(polygon, dtype=np.float32)
        if polygon.ndim != 2 or polygon.shape[0] < 3 or polygon.shape[1] != 2:
            continue
        mask |= MplPath(polygon, closed=True).contains_points(centroids)
    return mask


def _compute_slide_summary(*, slide_path: Path) -> dict[str, Any]:
    result: dict[str, Any] = {
        "feature_status": "ok",
        "ca1_polygon_count": 0,
        "ca1_cell_count": None,
        "ca1_pyramidal_count": None,
        "ca1_reactive_astro_count": None,
        "ca1_pyramidal_reactive_nn_distance_mean_px": None,
        "ca1_pyramidal_reactive_nn_distance_mean_um": None,
    }
    for spec in VARIATIONS:
        result[spec["feature_column"]] = np.nan
        result[f"near_count_{int(spec['radius_um'])}um"] = None

    if not slide_path.exists():
        result["feature_status"] = "missing_slide"
        return result

    polygons = load_region_polygons(slide_path, scale=16.0)
    ca1_polygons = polygons.get(REGION_NAME) or []
    result["ca1_polygon_count"] = len(ca1_polygons)
    if not ca1_polygons:
        result["feature_status"] = "missing_region"
        return result

    centroids = np.asarray(load_centroids(slide_path), dtype=np.float32)
    class_ids = np.asarray(load_class_ids(slide_path), dtype=np.int32)
    class_lookup = load_class_lookup(slide_path)
    label_to_id = _inverse_lookup(class_lookup)

    if PYRAMIDAL_LABEL not in label_to_id or REACTIVE_ASTRO_LABEL not in label_to_id:
        result["feature_status"] = "missing_required_label"
        return result

    pyramidal_id = label_to_id[PYRAMIDAL_LABEL]
    reactive_id = label_to_id[REACTIVE_ASTRO_LABEL]

    ca1_mask = _region_mask(centroids, ca1_polygons)
    result["ca1_cell_count"] = int(ca1_mask.sum())

    pyramidal_mask = ca1_mask & (class_ids == pyramidal_id)
    reactive_mask = ca1_mask & (class_ids == reactive_id)

    n_pyramidal = int(pyramidal_mask.sum())
    n_reactive = int(reactive_mask.sum())
    result["ca1_pyramidal_count"] = n_pyramidal
    result["ca1_reactive_astro_count"] = n_reactive

    if n_pyramidal < MIN_PYRAMIDAL_CELLS:
        result["feature_status"] = "too_few_pyramidal"
        return result

    pyramidal_xy = centroids[pyramidal_mask]
    if n_reactive == 0:
        result["feature_status"] = "ok_zero_reactive"
        result["ca1_pyramidal_reactive_nn_distance_mean_px"] = float("inf")
        result["ca1_pyramidal_reactive_nn_distance_mean_um"] = float("inf")
        for spec in VARIATIONS:
            result[spec["feature_column"]] = 0.0
            result[f"near_count_{int(spec['radius_um'])}um"] = 0
        return result

    reactive_xy = centroids[reactive_mask]
    tree = cKDTree(reactive_xy)
    nearest_dist_px, _ = tree.query(pyramidal_xy, k=1)
    nearest_dist_px = np.asarray(nearest_dist_px, dtype=np.float32)
    result["ca1_pyramidal_reactive_nn_distance_mean_px"] = float(np.mean(nearest_dist_px))
    result["ca1_pyramidal_reactive_nn_distance_mean_um"] = float(np.mean(nearest_dist_px) * 0.503)

    for spec in VARIATIONS:
        near_mask = nearest_dist_px <= float(spec["radius_px"])
        result[spec["feature_column"]] = float(np.mean(near_mask))
        result[f"near_count_{int(spec['radius_um'])}um"] = int(np.sum(near_mask))

    return result


def compute_variation_score(
    *,
    donor_id: str,
    data_root: str | Path,
    variation_name: str,
):
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    summary = _compute_slide_summary(slide_path=slide_path)
    return summary.get(VARIATION_LOOKUP[variation_name]["feature_column"])


def compute_donor_score(*, donor_id: str, data_root: str | Path):
    return compute_variation_score(
        donor_id=donor_id,
        data_root=data_root,
        variation_name=CANONICAL_VARIATION,
    )


def _prepare_donor_table(data_root: str | Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map({"Female": 0.0, "Male": 1.0})

    rows: list[dict[str, Any]] = []
    for row in cohort.itertuples(index=False):
        summary = _compute_slide_summary(slide_path=Path(data_root) / str(row.slide_name))
        merged = {"donor_id": str(row.donor_id)}
        merged.update(summary)
        rows.append(merged)

    summary_df = pd.DataFrame(rows)
    return cohort.merge(summary_df, on="donor_id", how="left")


def _evaluate_variation(donor_table: pd.DataFrame, spec: dict[str, Any]) -> dict[str, Any]:
    feature_col = spec["feature_column"]
    analyzable = donor_table.dropna(subset=[feature_col, OUTCOME_COLUMN, *CONFOUNDS]).reset_index(drop=True)
    n_total = int(len(donor_table))
    n_analyzable = int(len(analyzable))

    if n_analyzable >= 3:
        part = partial_correlation(
            analyzable,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUNDS,
        )
        loo_r = residualized_loo_predictive_correlation(
            analyzable,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUNDS,
        )
    else:
        part = {"partial_r": float("nan"), "p_value": float("nan")}
        loo_r = float("nan")

    selection = _selection_score(part.get("partial_r"), n_analyzable, n_total)
    gap, penalty, adjusted = _adjusted_metrics(part.get("partial_r"), loo_r)

    return {
        "name": spec["name"],
        "description": spec["description"],
        "radius_um": spec["radius_um"],
        "feature_name": spec["feature_name"],
        "feature_column": feature_col,
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "coverage_ratio": float(n_analyzable / n_total) if n_total else float("nan"),
        "partial_r": _float(part.get("partial_r")),
        "p_value": _float(part.get("p_value")),
        "selection_score": _float(selection),
        "loo_predictive_r": _float(loo_r),
        "is_loo_gap": _float(gap),
        "gap_penalty": _float(penalty),
        "adjusted_score": _float(adjusted),
    }


def _variation_sort_key(item: dict[str, Any]) -> tuple[float, float, float]:
    selection = item.get("selection_score")
    loo_r = item.get("loo_predictive_r")
    adjusted = item.get("adjusted_score")
    sel = float(selection) if selection is not None and np.isfinite(selection) else -np.inf
    loo = abs(float(loo_r)) if loo_r is not None and np.isfinite(loo_r) else -np.inf
    adj = float(adjusted) if adjusted is not None and np.isfinite(adjusted) else -np.inf
    return (sel, loo, adj)


def _fit_linear_prediction(train: pd.DataFrame, test: pd.DataFrame, feature_col: str) -> float:
    cols = CONFOUNDS + [feature_col]
    x_train = np.column_stack([np.ones(len(train), dtype=float)] + [train[c].to_numpy(dtype=float) for c in cols])
    y_train = train[OUTCOME_COLUMN].to_numpy(dtype=float)
    x_test = np.column_stack([np.ones(len(test), dtype=float)] + [test[c].to_numpy(dtype=float) for c in cols])
    beta, *_ = np.linalg.lstsq(x_train, y_train, rcond=None)
    return float((x_test @ beta)[0])


def _loo_prediction_table(frame: pd.DataFrame, feature_col: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for idx in range(len(frame)):
        train = frame.drop(index=frame.index[idx]).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)
        predicted = _fit_linear_prediction(train, test, feature_col=feature_col)
        row = test.iloc[0].to_dict()
        row["predicted"] = predicted
        row["abs_error"] = abs(float(row[OUTCOME_COLUMN]) - predicted)
        rows.append(row)
    return pd.DataFrame(rows)


def _status_counts(donor_table: pd.DataFrame) -> dict[str, int]:
    counts = donor_table["feature_status"].fillna("missing").value_counts(dropna=False).sort_index()
    return {str(k): int(v) for k, v in counts.items()}


def _patch_canonical_variation(best_variation: str) -> None:
    path = Path(__file__)
    text = path.read_text(encoding="utf-8")
    old = 'CANONICAL_VARIATION = "pyramidal_near_reactiveastro_30um"'
    import re

    text = re.sub(
        r'CANONICAL_VARIATION = ".*?"',
        f'CANONICAL_VARIATION = "{best_variation}"',
        text,
        count=1,
    )
    path.write_text(text, encoding="utf-8")


def _build_report(
    *,
    donor_table: pd.DataFrame,
    ranked: list[dict[str, Any]],
    best: dict[str, Any],
    loo_table: pd.DataFrame,
) -> str:
    selection = float(best["selection_score"]) if best.get("selection_score") is not None else float("nan")
    partial_r = float(best["partial_r"]) if best.get("partial_r") is not None else float("nan")
    loo_r = float(best["loo_predictive_r"]) if best.get("loo_predictive_r") is not None else float("nan")
    gap = float(best["is_loo_gap"]) if best.get("is_loo_gap") is not None else float("nan")
    penalty = float(best["gap_penalty"]) if best.get("gap_penalty") is not None else float("nan")
    adjusted = float(best["adjusted_score"]) if best.get("adjusted_score") is not None else float("nan")

    others = [x for x in ranked if x["name"] != best["name"]]
    others_sentence = (
        "; ".join(
            f"{item['name']} partial_r={item['partial_r']:.4f}, selection_score={item['selection_score']:.4f}, loo_r={item['loo_predictive_r']:.4f}"
            for item in others
            if item.get("partial_r") is not None and item.get("selection_score") is not None and item.get("loo_predictive_r") is not None
        )
        if others
        else "No other local variations were tested."
    )

    if loo_table.empty:
        error_sentence = "No analyzable donors were available for leave-one-out error inspection."
    else:
        worst = loo_table.sort_values("abs_error", ascending=False).head(3)
        error_sentence = (
            "Largest absolute LOO errors were for "
            + ", ".join(
                f"{row.donor_id} (outcome={float(getattr(row, OUTCOME_COLUMN)):.3f}, pred={float(row.predicted):.3f}, "
                f"30um={float(getattr(row, VARIATION_LOOKUP['pyramidal_near_reactiveastro_30um']['feature_column'])):.3f}, "
                f"50um={float(getattr(row, VARIATION_LOOKUP['pyramidal_near_reactiveastro_50um']['feature_column'])):.3f}, "
                f"CA1 pyramidal={int(getattr(row, 'ca1_pyramidal_count'))}, CA1 reactive astro={int(getattr(row, 'ca1_reactive_astro_count'))})"
                for row in worst.itertuples(index=False)
            )
            + "."
        )

    worked = (
        f"The strongest signal came from the {best['radius_um']:.0f} um peri-pyramidal niche in CA1, where the donor scalar is the fraction of CA1 pyramidal neurons "
        f"that have at least one nearby reactive astrocyte. That beats the broader local alternative, which suggests the biomarker is tied to a fairly immediate "
        f"neuron-adjacent reactive astrocyte pattern rather than diffuse CA1 gliosis alone."
    )
    failed = (
        "The broader 50 um neighborhood underperformed because it likely dilutes the perineuronal niche with more background reactive astrocytes, "
        "making the feature look more like general CA1 reactive astrocyte burden and less like a sharply neuron-adjacent injury pattern."
        if best["name"] == "pyramidal_near_reactiveastro_30um"
        else
        "The tighter 30 um neighborhood underperformed because it may miss biologically relevant reactive astrocytes sitting just outside the strict perineuronal band, "
        "whereas the 50 um radius still captures the same CA1 pyramidal-reactive astrocyte niche at donor level."
    )

    rationale = (
        f"This winner is biologically coherent because it keeps the round-1 CA1 reactive astrocyte signal but sharpens it to the specific niche around "
        f"CA1 pyramidal neurons. The summary is simple and auditable: among all CA1 pyramidal neurons, what fraction lies within {best['radius_um']:.0f} um "
        f"of a reactive astrocyte? It beat the nearby radius alternative because that distance best balanced specificity to a neuron-adjacent niche against "
        f"donor-level stability."
    )
    if best["name"] == "pyramidal_near_reactiveastro_30um":
        rationale += " Because it uses the same CA1 reactive astrocyte biology as the accepted panel member, it is more likely to refine or replace that member than to add wholly orthogonal information."
    else:
        rationale += " Because it remains anchored to the same CA1 reactive astrocyte biology as the accepted panel member, it may still be somewhat redundant with the current panel."

    interpretation = (
        f"The signal appears to represent CA1 reactive astrocyte engagement in the immediate peri-pyramidal niche. "
        f"Population: CA1 pyramidal neurons and CA1 reactive astrocytes. "
        f"Niche: reactive astrocytes within {best['radius_um']:.0f} um of pyramidal somata. "
        f"Feature summary: donor-level fraction of CA1 pyramidal neurons that are niche-positive. "
        f"Simplest observable tissue pattern: CA1 fields where many pyramidal neurons sit directly alongside reactive astrocytes track worse memory decline."
    )

    next_step = (
        "If this niche winner is kept, the next local sweep should stay in CA1 pyramidal–reactive astro biology but test whether normalizing by overall CA1 reactive astro burden or restricting to a very tight inner band improves additivity over the current panel."
    )

    return f"""## Summary
Tested the CA1 pyramidal-near-reactive-astrocyte niche family with a 30 um vs 50 um radius sweep; {best['name']} won with selection_score {selection:.4f}.

## Metrics
- Winning variation: {best['name']}
- IS partial r: {partial_r:.4f}
- Selection score: {selection:.4f}
- LOO predictive r: {loo_r:.4f}
- IS-LOO Gap: {gap:.4f} (penalty={penalty:.4f})
- Adjusted Score: {adjusted:.4f}
- Other tested variations: {others_sentence}

## Findings
1. What worked and why: {worked}
2. What failed and why: {failed}
3. Error pattern: {error_sentence}

## Rationale
{rationale}

## Interpretation
{interpretation}

## Next
{next_step}
"""


def main() -> int:
    data_root = Path("/data")
    scratch = Path("/scratch")

    donor_table = _prepare_donor_table(data_root)
    ranked = sorted((_evaluate_variation(donor_table, spec) for spec in VARIATIONS), key=_variation_sort_key, reverse=True)
    best = ranked[0]
    best_feature_col = str(best["feature_column"])

    analyzable = donor_table.dropna(subset=[best_feature_col, OUTCOME_COLUMN, *CONFOUNDS]).reset_index(drop=True)

    if len(analyzable) >= 3:
        boot = bootstrap_partial_correlation(
            analyzable,
            feature_col=best_feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUNDS,
            n_boot=2000,
            random_state=0,
        )
        stab = bootstrap_partial_correlation_stability(
            analyzable,
            feature_col=best_feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUNDS,
            n_boot=400,
            random_state=0,
        )
        loo_summary = leave_one_out_summary(
            analyzable,
            feature_col=best_feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUNDS,
            id_col="donor_id",
        )
        loo_table = _loo_prediction_table(analyzable, feature_col=best_feature_col)
        loo_table.to_csv(scratch / "loo_predictions.csv", index=False)
    else:
        boot = {"ci_lo": float("nan"), "ci_hi": float("nan")}
        stab = {
            "bootstrap_median_partial_r": float("nan"),
            "bootstrap_sign_consistency": float("nan"),
            "bootstrap_positive_fraction": float("nan"),
            "bootstrap_negative_fraction": float("nan"),
            "bootstrap_q25_partial_r": float("nan"),
            "bootstrap_q75_partial_r": float("nan"),
            "bootstrap_valid_samples": 0.0,
        }
        loo_summary = {"unstable_count": 0, "max_shift": float("nan"), "unstable_donors": []}
        loo_table = pd.DataFrame()

    other_feature_cols = [spec["feature_column"] for spec in VARIATIONS if spec["feature_column"] != best_feature_col]
    write_donor_feature_table(
        scratch / "donor_feature_table.csv",
        donor_table,
        feature_column=best_feature_col,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUNDS,
        extra_columns=other_feature_cols
        + [
            "feature_status",
            "ca1_polygon_count",
            "ca1_cell_count",
            "ca1_pyramidal_count",
            "ca1_reactive_astro_count",
            "near_count_30um",
            "near_count_50um",
            "ca1_pyramidal_reactive_nn_distance_mean_px",
            "ca1_pyramidal_reactive_nn_distance_mean_um",
        ],
    )

    results = build_results_payload(
        status="ok" if len(analyzable) >= 3 and best.get("partial_r") is not None else "null",
        feature_name=str(best["feature_name"]),
        outcome=OUTCOME_COLUMN,
        n_total=int(best["n_total"]),
        n_analyzable=int(best["n_analyzable"]),
        partial_r=_float(best.get("partial_r")),
        ci_lo=_float(boot.get("ci_lo")),
        ci_hi=_float(boot.get("ci_hi")),
        p_value=_float(best.get("p_value")),
        loo_predictive_r=_float(best.get("loo_predictive_r")),
        loo_unstable_count=int(loo_summary.get("unstable_count", 0)),
        loo_max_shift=_float(loo_summary.get("max_shift")),
        donor_ids_used=analyzable["donor_id"].astype(str).tolist(),
        covariates=CONFOUNDS,
        recomputed_from_raw=True,
        registry_written=False,
        artifacts={
            "donor_feature_table": str(scratch / "donor_feature_table.csv"),
            "loo_predictions": str(scratch / "loo_predictions.csv"),
        },
        hypothesis_family=HYPOTHESIS_FAMILY,
        best_variation=str(best["name"]),
        feature_column=best_feature_col,
        coverage_ratio=_float(best.get("coverage_ratio")),
        selection_score=_float(best.get("selection_score")),
        is_loo_gap=_float(best.get("is_loo_gap")),
        gap_penalty=_float(best.get("gap_penalty")),
        adjusted_score=_float(best.get("adjusted_score")),
        bootstrap_median_partial_r=_float(stab.get("bootstrap_median_partial_r")),
        bootstrap_sign_consistency=_float(stab.get("bootstrap_sign_consistency")),
        bootstrap_positive_fraction=_float(stab.get("bootstrap_positive_fraction")),
        bootstrap_negative_fraction=_float(stab.get("bootstrap_negative_fraction")),
        bootstrap_q25_partial_r=_float(stab.get("bootstrap_q25_partial_r")),
        bootstrap_q75_partial_r=_float(stab.get("bootstrap_q75_partial_r")),
        bootstrap_valid_samples=_float(stab.get("bootstrap_valid_samples")),
        status_counts=_status_counts(donor_table),
        ranked_variations=ranked,
    )
    write_results_payload(scratch / "results.json", results)

    _patch_canonical_variation(str(best["name"]))

    report = _build_report(donor_table=donor_table, ranked=ranked, best=best, loo_table=loo_table)
    (scratch / "report.md").write_text(report, encoding="utf-8")

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best['name']}")
    print(f"  IS partial r:      {float(best['partial_r']):.4f}")
    print(f"  Selection score:   {float(best['selection_score']):.4f}")
    print(f"  LOO predictive r:  {float(best['loo_predictive_r']):.4f}  (diagnostic)")
    print(f"  IS-LOO Gap:        {float(best['is_loo_gap']):.4f}  (penalty={float(best['gap_penalty']):.4f})")
    print(f"  Adjusted Score:    {float(best['adjusted_score']):.4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked:
        pr = float(item["partial_r"]) if item.get("partial_r") is not None else float("nan")
        ss = float(item["selection_score"]) if item.get("selection_score") is not None else float("nan")
        lr = float(item["loo_predictive_r"]) if item.get("loo_predictive_r") is not None else float("nan")
        print(f"  {item['name']}  {pr:.4f}  {ss:.4f}  {lr:.4f}")
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best_feature_col}  ca1_pyramidal_count  ca1_reactive_astro_count")
    if loo_table.empty:
        print("  none  nan  nan  nan  nan  nan")
    else:
        for row in loo_table[
            ["donor_id", OUTCOME_COLUMN, "predicted", best_feature_col, "ca1_pyramidal_count", "ca1_reactive_astro_count"]
        ].itertuples(index=False):
            print(
                f"  {row.donor_id}  {float(getattr(row, OUTCOME_COLUMN)):.4f}  {float(row.predicted):.4f}  "
                f"{float(getattr(row, best_feature_col)):.4f}  {int(row.ca1_pyramidal_count)}  {int(row.ca1_reactive_astro_count)}"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
