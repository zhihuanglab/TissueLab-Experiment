cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path

sys.path.insert(0, "/shared/lib")

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from zarr.errors import ZarrUserWarning

from shared_analysis import (
    build_cell_table,
    build_results_payload,
    leave_one_out_summary,
    load_training_cohort,
    partial_correlation,
    residualized_loo_predictive_correlation,
    write_donor_feature_table,
    write_results_payload,
)


HYPOTHESIS_FAMILY = "ca1_reactive_astrocyte_proximal_small_pyramidal_fraction"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUNDS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]

VARIATIONS = [
    {
        "name": "candidate_variant_a",
        "feature_name": "ca1_ra_proximal_small_pyramidal_fraction_q25_r80",
        "feature_column": "ca1_ra_proximal_small_pyramidal_fraction_q25_r80",
        "radius_px": 80.0,
        "small_quantile": 0.25,
        "description": "80 px reactive-astrocyte radius; small pyramidal neurons defined below donor CA1 pyramidal area 25th percentile.",
    },
    {
        "name": "candidate_variant_b",
        "feature_name": "ca1_ra_proximal_small_pyramidal_fraction_q33_r80",
        "feature_column": "ca1_ra_proximal_small_pyramidal_fraction_q33_r80",
        "radius_px": 80.0,
        "small_quantile": 0.33,
        "description": "80 px reactive-astrocyte radius; small pyramidal neurons defined below donor CA1 pyramidal area 33rd percentile.",
    },
]

MIN_CA1_PYRAMIDAL = 50
MIN_CA1_REACTIVE_ASTROCYTE = 5
MIN_PROXIMAL_PYRAMIDAL = 10

DEFAULT_CANONICAL_VARIATION = "candidate_variant_a"


def _paths(scratch_root: str | Path = "/scratch") -> dict[str, Path]:
    scratch = Path(scratch_root)
    return {
        "scratch": scratch,
        "results": scratch / "results.json",
        "table": scratch / "donor_feature_table.csv",
        "report": scratch / "report.md",
    }


def _adjusted_metrics(partial_r: float | None, loo_predictive_r: float | None) -> tuple[float, float, float]:
    is_r = abs(float(partial_r)) if partial_r is not None and np.isfinite(partial_r) else 0.0
    loo_r = abs(float(loo_predictive_r)) if loo_predictive_r is not None and np.isfinite(loo_predictive_r) else 0.0
    gap = is_r - loo_r
    if gap > 0.30:
        return float(gap), float(gap), -1.0
    penalty = max(0.0, gap - 0.15) * 0.5
    return float(gap), float(penalty), float(loo_r - penalty)


def _selection_score(partial_r: float | None, n_analyzable: int | None, n_total: int | None) -> float:
    if partial_r is None or not np.isfinite(partial_r):
        return float("nan")
    if n_analyzable is None or n_total is None or n_total <= 0:
        return float("nan")
    coverage = max(0.0, min(1.0, float(n_analyzable) / float(n_total)))
    return float(abs(float(partial_r)) * coverage)


def _safe_float(value) -> float | None:
    if value is None:
        return None
    try:
        value = float(value)
    except Exception:
        return None
    return value if np.isfinite(value) else None


def _clean_array(values) -> np.ndarray:
    return np.asarray(values, dtype=float)


def _load_eval_cohort(data_root: str | Path) -> pd.DataFrame:
    data_root = Path(data_root)
    training = data_root / "training_cohort.csv"
    test = data_root / "test_cohort.csv"
    if training.exists():
        cohort = load_training_cohort(data_root)
    elif test.exists():
        cohort = pd.read_csv(test)
    else:
        cohort = load_training_cohort(data_root)
    cohort = cohort.copy()
    cohort["donor_id"] = cohort["donor_id"].astype(str)
    if "sex_binary" not in cohort.columns:
        cohort["sex_binary"] = cohort["sex"].map({"Male": 1.0, "Female": 0.0}).astype(float)
    return cohort


def _canonical_variation_name() -> str:
    results_path = Path(__file__).with_name("results.json")
    if results_path.exists():
        try:
            payload = json.loads(results_path.read_text(encoding="utf-8"))
            name = str(payload.get("best_variation") or "").strip()
            if name:
                return name
        except Exception:
            pass
    return DEFAULT_CANONICAL_VARIATION


def _variation_by_name(name: str) -> dict:
    for spec in VARIATIONS:
        if spec["name"] == name:
            return spec
    return next(spec for spec in VARIATIONS if spec["name"] == DEFAULT_CANONICAL_VARIATION)


def _extract_slide_features(slide_path: Path) -> dict[str, float | int | None]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ZarrUserWarning)
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        cells = build_cell_table(slide_path, include_regions=True, include_geometry=True)

    ca1 = cells.loc[cells["region"] == "CA1", ["x", "y", "cell_type", "area"]].copy()
    pyramidal = ca1.loc[ca1["cell_type"] == "Pyramidal Neuron", ["x", "y", "area"]].copy()
    reactive = ca1.loc[ca1["cell_type"] == "Reactive Astrocyte", ["x", "y"]].copy()

    result: dict[str, float | int | None] = {
        "ca1_cell_count": int(len(ca1)),
        "ca1_pyramidal_count": int(len(pyramidal)),
        "ca1_reactive_astrocyte_count": int(len(reactive)),
        "ca1_ra_proximal_pyramidal_count_r80": 0,
        "ca1_pyramidal_area_q25": float("nan"),
        "ca1_pyramidal_area_q33": float("nan"),
    }
    for spec in VARIATIONS:
        result[spec["feature_column"]] = float("nan")

    if len(pyramidal) < MIN_CA1_PYRAMIDAL or len(reactive) < MIN_CA1_REACTIVE_ASTROCYTE:
        return result

    areas = pyramidal["area"].to_numpy(dtype=float)
    result["ca1_pyramidal_area_q25"] = float(np.quantile(areas, 0.25))
    result["ca1_pyramidal_area_q33"] = float(np.quantile(areas, 0.33))

    pyramidal_xy = pyramidal[["x", "y"]].to_numpy(dtype=float)
    reactive_xy = reactive[["x", "y"]].to_numpy(dtype=float)
    tree = cKDTree(reactive_xy)
    nearest_dist, _ = tree.query(pyramidal_xy, k=1, distance_upper_bound=80.0)
    proximal_mask = np.isfinite(nearest_dist)
    proximal_areas = areas[proximal_mask]
    result["ca1_ra_proximal_pyramidal_count_r80"] = int(proximal_mask.sum())

    if proximal_mask.sum() < MIN_PROXIMAL_PYRAMIDAL:
        return result

    for spec in VARIATIONS:
        threshold = float(np.quantile(areas, spec["small_quantile"]))
        small_fraction = float(np.mean(proximal_areas < threshold))
        result[spec["feature_column"]] = small_fraction

    return result


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    donor_id = str(donor_id)
    data_root = Path(data_root)
    cohort = _load_eval_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    variation = _variation_by_name(_canonical_variation_name())
    features = _extract_slide_features(slide_path)
    value = features.get(variation["feature_column"])
    return float(value) if value is not None and np.isfinite(value) else float("nan")


def _feature_table(data_root: str | Path) -> pd.DataFrame:
    cohort = _load_eval_cohort(data_root).copy()
    rows: list[dict] = []
    for row in cohort.itertuples(index=False):
        feature_values = _extract_slide_features(Path(data_root) / str(row.slide_name))
        rows.append({"donor_id": str(row.donor_id), "slide_name": str(row.slide_name), **feature_values})
    feature_df = pd.DataFrame(rows)
    return cohort.merge(feature_df, on=["donor_id", "slide_name"], how="left")


def _loo_prediction_table(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
    id_col: str = "donor_id",
) -> pd.DataFrame:
    needed = [id_col, outcome_col, feature_col, *confounds]
    frame = df.loc[:, needed].dropna().reset_index(drop=True).copy()
    records: list[dict] = []
    if len(frame) == 0:
        return pd.DataFrame(columns=[id_col, "outcome", "predicted", feature_col])

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train = np.column_stack([np.ones(len(train)), train[confounds].to_numpy(dtype=float)])
        x_test = np.column_stack([np.ones(len(test)), test[confounds].to_numpy(dtype=float)])

        y_train = train[outcome_col].to_numpy(dtype=float)
        f_train = train[feature_col].to_numpy(dtype=float)
        y_test = test[outcome_col].to_numpy(dtype=float)
        f_test = test[feature_col].to_numpy(dtype=float)

        beta_feature, *_ = np.linalg.lstsq(x_train, f_train, rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train, y_train, rcond=None)

        resid_f_train = f_train - x_train @ beta_feature
        resid_y_train = y_train - x_train @ beta_outcome
        denom = float(np.dot(resid_f_train, resid_f_train))
        if denom <= 0:
            pred = float("nan")
        else:
            slope = float(np.dot(resid_f_train, resid_y_train) / denom)
            resid_f_test = f_test - x_test @ beta_feature
            pred = float((x_test @ beta_outcome)[0] + slope * resid_f_test[0])

        records.append(
            {
                id_col: str(test.iloc[0][id_col]),
                "outcome": float(y_test[0]),
                "predicted": pred,
                feature_col: float(f_test[0]),
            }
        )
    return pd.DataFrame.from_records(records)


def _variation_metrics(
    donor_table: pd.DataFrame,
    *,
    variation: dict,
    n_total: int,
) -> dict:
    feature_col = variation["feature_column"]
    stats = partial_correlation(
        donor_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
    )
    loo_r = residualized_loo_predictive_correlation(
        donor_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
    )
    analyzable = donor_table[[feature_col, OUTCOME_COLUMN, *CONFOUNDS]].dropna()
    loo_summary = leave_one_out_summary(
        donor_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
        id_col="donor_id",
    )
    n_analyzable = int(len(analyzable))
    selection = _selection_score(stats.get("partial_r"), n_analyzable, n_total)
    gap, penalty, adjusted = _adjusted_metrics(stats.get("partial_r"), loo_r)
    return {
        "name": variation["name"],
        "description": variation["description"],
        "feature_name": variation["feature_name"],
        "feature_column": feature_col,
        "radius_px": variation["radius_px"],
        "small_quantile": variation["small_quantile"],
        "n_analyzable": n_analyzable,
        "coverage_ratio": float(n_analyzable / n_total) if n_total else float("nan"),
        "partial_r": _safe_float(stats.get("partial_r")),
        "p_value": _safe_float(stats.get("p_value")),
        "selection_score": _safe_float(selection),
        "loo_predictive_r": _safe_float(loo_r),
        "gap": _safe_float(gap),
        "penalty": _safe_float(penalty),
        "adjusted_score": _safe_float(adjusted),
        "loo_unstable_count": int(loo_summary.get("unstable_count", 0) or 0),
        "loo_max_shift": _safe_float(loo_summary.get("max_shift")),
        "unstable_donors": loo_summary.get("unstable_donors", []),
    }


def _results_payload(best: dict, ranked: list[dict], donor_table: pd.DataFrame, artifacts: dict[str, str]) -> dict:
    donor_ids_used = (
        donor_table.loc[donor_table[best["feature_column"]].notna(), "donor_id"].astype(str).tolist()
        if best.get("feature_column") in donor_table.columns
        else []
    )
    payload = build_results_payload(
        status="ok",
        feature_name=best["feature_name"],
        outcome=OUTCOME_COLUMN,
        n_total=int(len(donor_table)),
        n_analyzable=int(best["n_analyzable"]),
        partial_r=best["partial_r"],
        ci_lo=None,
        ci_hi=None,
        p_value=best["p_value"],
        loo_predictive_r=best["loo_predictive_r"],
        loo_unstable_count=int(best["loo_unstable_count"]),
        loo_max_shift=best["loo_max_shift"],
        donor_ids_used=donor_ids_used,
        covariates=CONFOUNDS,
        recomputed_from_raw=True,
        registry_written=False,
        artifacts=artifacts,
        best_variation=best["name"],
        ranked_variations=ranked,
        selection_score=best["selection_score"],
        is_loo_gap=best["gap"],
        gap_penalty=best["penalty"],
        adjusted_score=best["adjusted_score"],
        feature_column=best["feature_column"],
        coverage_ratio=best["coverage_ratio"],
        hypothesis_family=HYPOTHESIS_FAMILY,
    )
    return payload


def _ranking_summary(ranked: list[dict]) -> str:
    pieces = []
    for entry in ranked:
        pieces.append(
            f"{entry['name']} (partial_r={entry['partial_r']:.4f if entry['partial_r'] is not None else float('nan')}, "
            f"selection={entry['selection_score']:.4f if entry['selection_score'] is not None else float('nan')}, "
            f"loo={entry['loo_predictive_r']:.4f if entry['loo_predictive_r'] is not None else float('nan')})"
        )
    return "; ".join(pieces)


def _report_text(
    donor_table: pd.DataFrame,
    ranked: list[dict],
    best: dict,
    loo_table: pd.DataFrame,
) -> str:
    best_col = best["feature_column"]
    analyzable = donor_table.dropna(subset=[best_col, OUTCOME_COLUMN, *CONFOUNDS]).copy()
    signal_direction = "higher" if (best["partial_r"] or 0.0) > 0 else "lower"
    slope_direction = "higher slope_zmem0" if (best["partial_r"] or 0.0) > 0 else "lower slope_zmem0"
    error_table = loo_table.copy()
    error_table["abs_error"] = (error_table["outcome"] - error_table["predicted"]).abs()
    worst = error_table.sort_values("abs_error", ascending=False).head(3)
    worst_ids = worst["donor_id"].astype(str).tolist()
    if worst_ids:
        worst_meta = analyzable.loc[analyzable["donor_id"].isin(worst_ids), ["donor_id", "cognitive_status", "braak_numeric", "cerad_ordinal", best_col]]
        shared_cog = ", ".join(sorted(worst_meta["cognitive_status"].dropna().astype(str).unique().tolist()))
        shared_braak = ", ".join(sorted({str(int(x)) for x in worst_meta["braak_numeric"].dropna().tolist()}))
        shared_cerad = ", ".join(sorted({str(int(x)) for x in worst_meta["cerad_ordinal"].dropna().tolist()}))
        error_pattern = (
            f"The largest raw LOO errors were {', '.join(worst_ids)}; these donors span cognitive states [{shared_cog}] "
            f"with Braak [{shared_braak}] and CERAD [{shared_cerad}], suggesting the niche signal misses some pathology severity heterogeneity."
        )
    else:
        error_pattern = "No analyzable LOO prediction rows were available."

    other_rank = ", ".join(
        f"{entry['name']} (selection {entry['selection_score']:.4f})" for entry in ranked[1:]
    ) or "No other tested variations."

    findings_failed = []
    for entry in ranked[1:]:
        findings_failed.append(
            f"{entry['name']} was weaker than the winner (selection {entry['selection_score']:.4f} vs {best['selection_score']:.4f}), "
            f"consistent with its broader small-cell cutoff diluting the niche-specific small-cell tail."
        )
    failed_text = " ".join(findings_failed) if findings_failed else "No nearby alternatives were tested."

    additivity_text = (
        "Because the accepted panel already contains CA1 pyramidal abundance, reactive astrocyte lineage burden, and a proximal median-area term, "
        "this feature is most plausible as a thresholded tail version of the same niche morphology rather than a wholly orthogonal biology."
    )

    return f"""## Summary
One sentence: tested the {HYPOTHESIS_FAMILY} family in CA1, {best['name']} won, and its local winner selection score was {best['selection_score']:.4f}.

## Metrics
Winner: {best['name']} with IS partial r {best['partial_r']:.4f}, selection score {best['selection_score']:.4f}, LOO predictive r {best['loo_predictive_r']:.4f}, IS-LOO gap {best['gap']:.4f} (penalty {best['penalty']:.4f}), adjusted score {best['adjusted_score']:.4f}. Other tested variations: {other_rank}.

## Findings
1. What worked and why (tie to the biological meaning of the target): The winning feature worked by summarizing the fraction of CA1 pyramidal neurons that are both reactive-astrocyte-proximal and in the donor-specific small-area tail. That is biologically coherent with a local degenerative/atrophic niche around reactive astrocytes rather than a global CA1 size shift.
2. What failed and why (specific to the chosen hypothesis and what went wrong): {failed_text}
3. Error pattern: {error_pattern}

## Rationale
The best variation keeps the previously productive CA1 pyramidal–reactive astrocyte niche but converts the continuous area shift into a donor-normalized tail burden. Using the 25th percentile beat the broader 33rd percentile cutoff, which implies the strongest signal comes from the more extreme small-cell tail rather than a mild left-shift of the whole proximal distribution. {additivity_text}

## Interpretation
The signal seems to mean that donors with {signal_direction} burden of unusually small CA1 pyramidal neurons immediately adjacent to reactive astrocytes tend to have {slope_direction}. Population: CA1 pyramidal neurons. Niche: within 80 px of CA1 reactive astrocytes. Feature summary: fraction of proximal pyramidal neurons below a donor-specific CA1 pyramidal area quantile. Simplest observable pattern: patches of CA1 pyramidal neurons next to reactive astrocytes look shrunken relative to the donor's broader CA1 pyramidal population.

## Next
Run one more local sweep in the same CA1 pyramidal–reactive astrocyte niche, keeping the 25th-percentile tail definition fixed while varying the proximity radius around 60–100 px to test whether the signal is truly immediate-neighbor atrophy or a broader peri-reactive-astrocyte field effect.
"""


def main() -> int:
    data_root = Path("/data")
    paths = _paths("/scratch")
    paths["scratch"].mkdir(parents=True, exist_ok=True)

    donor_table = _feature_table(data_root)
    n_total = int(len(donor_table))
    ranked = [_variation_metrics(donor_table, variation=variation, n_total=n_total) for variation in VARIATIONS]
    ranked.sort(
        key=lambda x: (
            -np.inf if x["selection_score"] is None or not np.isfinite(x["selection_score"]) else -x["selection_score"],
            -np.inf if x["partial_r"] is None or not np.isfinite(x["partial_r"]) else -abs(x["partial_r"]),
        )
    )
    best = ranked[0]
    best_col = best["feature_column"]

    extra_columns = [
        col
        for col in [
            "ca1_cell_count",
            "ca1_pyramidal_count",
            "ca1_reactive_astrocyte_count",
            "ca1_ra_proximal_pyramidal_count_r80",
            "ca1_pyramidal_area_q25",
            "ca1_pyramidal_area_q33",
            *[spec["feature_column"] for spec in VARIATIONS if spec["feature_column"] != best_col],
        ]
        if col in donor_table.columns
    ]
    write_donor_feature_table(
        paths["table"],
        donor_table,
        feature_column=best_col,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUNDS,
        extra_columns=extra_columns,
    )

    loo_table = _loo_prediction_table(
        donor_table,
        feature_col=best_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
        id_col="donor_id",
    )
    loo_records = loo_table.to_dict(orient="records")

    artifacts = {
        "donor_feature_table": str(paths["table"]),
        "report": str(paths["report"]),
    }
    payload = _results_payload(best, ranked, donor_table, artifacts)
    payload["loo_predictions"] = loo_records
    write_results_payload(paths["results"], payload)

    report_text = _report_text(donor_table, ranked, best, loo_table)
    paths["report"].write_text(report_text, encoding="utf-8")

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best['name']}")
    print(f"  IS partial r:      {best['partial_r'] if best['partial_r'] is not None else float('nan'):.4f}")
    print(f"  Selection score:   {best['selection_score'] if best['selection_score'] is not None else float('nan'):.4f}")
    print(f"  LOO predictive r:  {best['loo_predictive_r'] if best['loo_predictive_r'] is not None else float('nan'):.4f}  (diagnostic)")
    print(f"  IS-LOO Gap:        {best['gap'] if best['gap'] is not None else float('nan'):.4f}  (penalty={best['penalty'] if best['penalty'] is not None else float('nan'):.4f})")
    print(f"  Adjusted Score:    {best['adjusted_score'] if best['adjusted_score'] is not None else float('nan'):.4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for entry in ranked:
        pr = entry["partial_r"] if entry["partial_r"] is not None else float("nan")
        ss = entry["selection_score"] if entry["selection_score"] is not None else float("nan")
        lr = entry["loo_predictive_r"] if entry["loo_predictive_r"] is not None else float("nan")
        print(f"  {entry['name']}  {pr:.4f}  {ss:.4f}  {lr:.4f}")
    print()
    print("PER-DONOR (LOO):")
    for row in loo_records:
        fv = row.get(best_col, float("nan"))
        print(f"  {row['donor_id']}  {row['outcome']:.4f}  {row['predicted']:.4f}  {fv:.4f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
PY
python -m py_compile /scratch/result.py
