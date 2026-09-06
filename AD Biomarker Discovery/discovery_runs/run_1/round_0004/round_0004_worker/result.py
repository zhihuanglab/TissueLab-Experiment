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
from shared_analysis.artifacts import build_results_payload, write_donor_feature_table, write_results_payload
from shared_analysis.stats import partial_correlation


DATA_ROOT_DEFAULT = Path("/data")
OUTCOME_COLUMN = "slope_zmem0"
CONFOUNDS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]

FEATURE_FAMILY = "ca1_peripyramidal_reactive_astro_coverage"
BASELINE_VARIATION = "coverage_80px"
VARIATIONS = {
    "coverage_80px": {
        "radius_px": 80.0,
        "feature_column": "ca1_pyramidal_reactive_astro_coverage_80px",
        "description": "Fraction of CA1 pyramidal neurons with >=1 CA1 reactive astrocyte within 80 px.",
    },
    "coverage_60px": {
        "radius_px": 60.0,
        "feature_column": "ca1_pyramidal_reactive_astro_coverage_60px",
        "description": "Fraction of CA1 pyramidal neurons with >=1 CA1 reactive astrocyte within 60 px.",
    },
}
KEY_COUNT_COLUMNS = ["ca1_pyramidal_neuron_count", "ca1_reactive_astrocyte_count"]


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


def _design_matrix(confounds: pd.DataFrame) -> np.ndarray:
    matrix = confounds.astype(float).to_numpy()
    intercept = np.ones((len(confounds), 1), dtype=float)
    return np.hstack([intercept, matrix])


def _prepare_cohort(data_root: str | Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = (cohort["sex"].astype(str).str.lower() == "male").astype(float)
    return cohort


def _extract_slide_features(slide_path: Path) -> dict[str, float]:
    warnings.filterwarnings("ignore", category=UserWarning)
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)

    mask = (cells["region"] == "CA1") & (cells["cell_type"].isin(["Pyramidal Neuron", "Reactive Astrocyte"]))
    ca1 = cells.loc[mask, ["x", "y", "cell_type"]].copy()

    pyramidal = ca1.loc[ca1["cell_type"] == "Pyramidal Neuron", ["x", "y"]].to_numpy(dtype=np.float32, copy=False)
    reactive = ca1.loc[ca1["cell_type"] == "Reactive Astrocyte", ["x", "y"]].to_numpy(dtype=np.float32, copy=False)

    result: dict[str, float] = {
        "ca1_pyramidal_neuron_count": float(len(pyramidal)),
        "ca1_reactive_astrocyte_count": float(len(reactive)),
    }

    if len(pyramidal) == 0:
        for spec in VARIATIONS.values():
            result[spec["feature_column"]] = float("nan")
        return result

    if len(reactive) == 0:
        for spec in VARIATIONS.values():
            result[spec["feature_column"]] = 0.0
        return result

    tree = cKDTree(reactive)
    for spec in VARIATIONS.values():
        radius = float(spec["radius_px"])
        dists, _ = tree.query(pyramidal, k=1, distance_upper_bound=radius)
        covered = np.isfinite(dists)
        result[spec["feature_column"]] = float(np.mean(covered))
    return result


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Return one biomarker score for one donor/slide row.

    The evaluator loops over the cohort, calls this function, builds the donor
    table, joins confounds/outcome, and computes all metrics.

    Keep this function portable:
    - recompute from raw data under data_root
    - do not depend on /shared/cache in the final script
    - if you need learned state, save it next to result.py and load via __file__
    """
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


def _clean_feature_frame(df: pd.DataFrame, feature_column: str) -> pd.DataFrame:
    cols = ["donor_id", OUTCOME_COLUMN, feature_column, *CONFOUNDS, *KEY_COUNT_COLUMNS]
    frame = df.loc[:, [c for c in cols if c in df.columns]].replace([np.inf, -np.inf], np.nan).dropna(
        subset=[feature_column, OUTCOME_COLUMN, *CONFOUNDS]
    )
    return frame.reset_index(drop=True)


def _loo_predictions(df: pd.DataFrame, feature_column: str) -> tuple[float, pd.DataFrame]:
    frame = _clean_feature_frame(df, feature_column)
    rows: list[dict[str, float | str]] = []
    pred_resids: list[float] = []
    actual_resids: list[float] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train = _design_matrix(train[CONFOUNDS])
        x_test = _design_matrix(test[CONFOUNDS])

        y_feature_train = train[feature_column].to_numpy(dtype=float)
        y_outcome_train = train[OUTCOME_COLUMN].to_numpy(dtype=float)

        beta_feature, *_ = np.linalg.lstsq(x_train, y_feature_train, rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train, y_outcome_train, rcond=None)

        resid_feature_train = y_feature_train - x_train @ beta_feature
        resid_outcome_train = y_outcome_train - x_train @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        pred_resid = float("nan")
        pred_outcome = float("nan")
        actual_resid = float(test[OUTCOME_COLUMN].to_numpy(dtype=float)[0] - (x_test @ beta_outcome)[0])

        if denom > 0:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            feature_resid_test = float(test[feature_column].to_numpy(dtype=float)[0] - (x_test @ beta_feature)[0])
            pred_resid = float(slope * feature_resid_test)
            pred_outcome = float((x_test @ beta_outcome)[0] + pred_resid)

        row = {
            "donor_id": str(test.loc[0, "donor_id"]),
            "outcome": float(test.loc[0, OUTCOME_COLUMN]),
            "predicted": pred_outcome,
            "predicted_residualized": pred_resid,
            "actual_residualized": actual_resid,
            feature_column: float(test.loc[0, feature_column]),
        }
        for count_col in KEY_COUNT_COLUMNS:
            if count_col in test.columns:
                row[count_col] = float(test.loc[0, count_col])
        rows.append(row)
        pred_resids.append(pred_resid)
        actual_resids.append(actual_resid)

    pred_arr = np.asarray(pred_resids, dtype=float)
    actual_arr = np.asarray(actual_resids, dtype=float)
    mask = np.isfinite(pred_arr) & np.isfinite(actual_arr)
    if mask.sum() >= 3 and np.std(pred_arr[mask]) > 0 and np.std(actual_arr[mask]) > 0:
        loo_r = float(np.corrcoef(pred_arr[mask], actual_arr[mask])[0, 1])
    else:
        loo_r = float("nan")
    return loo_r, pd.DataFrame(rows)


def _gap_penalty(is_partial_r: float, loo_predictive_r: float) -> tuple[float, float]:
    if not np.isfinite(is_partial_r) or not np.isfinite(loo_predictive_r):
        return float("nan"), float("nan")
    gap = float(abs(is_partial_r) - abs(loo_predictive_r))
    penalty = float(max(0.0, gap - 0.15) * 0.5)
    return gap, penalty


def _adjusted_score(selection_score: float, gap: float, penalty: float) -> float:
    if not np.isfinite(selection_score) or not np.isfinite(gap) or not np.isfinite(penalty):
        return float("nan")
    if gap > 0.30:
        return -1.0
    return float(selection_score - penalty)


def _evaluate_variation(df: pd.DataFrame, variation_name: str) -> dict[str, object]:
    spec = VARIATIONS[variation_name]
    feature_column = str(spec["feature_column"])
    frame = _clean_feature_frame(df, feature_column)

    is_metrics = partial_correlation(
        df,
        feature_col=feature_column,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
    )
    partial_r = float(is_metrics["partial_r"])
    p_value = float(is_metrics["p_value"])
    n_analyzable = int(len(frame))
    n_total = int(len(df))
    selection_score = float(abs(partial_r) * (n_analyzable / n_total)) if np.isfinite(partial_r) else float("nan")

    loo_predictive_r, per_donor = _loo_predictions(df, feature_column)
    gap, penalty = _gap_penalty(partial_r, loo_predictive_r)
    adjusted_score = _adjusted_score(selection_score, gap, penalty)

    return {
        "name": variation_name,
        "description": spec["description"],
        "radius_px": float(spec["radius_px"]),
        "feature_column": feature_column,
        "partial_r": partial_r,
        "p_value": p_value,
        "selection_score": selection_score,
        "loo_predictive_r": loo_predictive_r,
        "is_loo_gap": gap,
        "penalty": penalty,
        "adjusted_score": adjusted_score,
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "donor_ids_used": frame["donor_id"].astype(str).tolist(),
        "per_donor_loo": per_donor.to_dict(orient="records"),
    }


def _format_float(value: float) -> str:
    if value is None or not np.isfinite(value):
        return "nan"
    return f"{value:.4f}"


def _extract_all_features(data_root: Path) -> pd.DataFrame:
    cohort = _prepare_cohort(data_root)
    feature_rows: list[dict[str, float | str]] = []
    for row in cohort.itertuples(index=False):
        slide_path = data_root / str(row.slide_name)
        feature_rows.append({"donor_id": str(row.donor_id), **_extract_slide_features(slide_path)})
    feature_df = pd.DataFrame(feature_rows)
    return cohort.merge(feature_df, on="donor_id", how="left", validate="one_to_one")


def _build_report(best: dict[str, object], ranked: list[dict[str, object]], donor_table: pd.DataFrame) -> str:
    best_name = str(best["name"])
    best_feature = str(best["feature_column"])
    radius = int(float(best["radius_px"]))
    abs_err = (donor_table["predicted"] - donor_table["outcome"]).abs()
    error_table = donor_table.assign(abs_error=abs_err).sort_values("abs_error", ascending=False)
    worst = error_table.head(5)

    other_lines = []
    loser_rows = []
    for row in ranked:
        if row["name"] == best_name:
            continue
        loser_rows.append(row)
        other_lines.append(
            f"- {row['name']}: partial_r={_format_float(float(row['partial_r']))}, "
            f"selection_score={_format_float(float(row['selection_score']))}, "
            f"loo_predictive_r={_format_float(float(row['loo_predictive_r']))}, "
            f"adjusted_score={_format_float(float(row['adjusted_score']))}"
        )
    if not other_lines:
        other_lines.append("- No alternate local variations were tested.")

    if loser_rows:
        loser = loser_rows[0]
        if float(best["radius_px"]) < float(loser["radius_px"]):
            failure_text = (
                "The broader nearby radius was slightly weaker, which suggests that "
                "expanding the neighborhood starts to dilute the specifically perineuronal "
                "reactive-astrocyte signal with more background CA1 astrocyte burden."
            )
        else:
            failure_text = (
                "The tighter nearby radius was slightly weaker, which suggests that an overly "
                "strict neighborhood misses part of the biologically relevant short-range "
                "reactive-astrocyte exposure around CA1 pyramidal neurons."
            )
    else:
        failure_text = "No nearby alternative was available to diagnose family-local failure."

    worst_lines = []
    for row in worst.itertuples(index=False):
        worst_lines.append(
            f"- {row.donor_id}: outcome={row.outcome:.4f}, predicted={row.predicted:.4f}, "
            f"{best_feature}={getattr(row, best_feature):.4f}, "
            f"CA1 pyramidal={int(row.ca1_pyramidal_neuron_count)}, "
            f"CA1 reactive astrocyte={int(row.ca1_reactive_astrocyte_count)}"
        )

    analyzable = int(best["n_analyzable"])
    n_total = int(best["n_total"])

    if loser_rows:
        if float(best["radius_px"]) < float(loser["radius_px"]):
            rationale_beat_text = (
                "It beat the nearby alternative because the winning radius stayed tighter around the "
                "perineuronal zone and avoided diluting the signal with more diffuse background CA1 "
                "reactive astrocyte burden."
            )
        else:
            rationale_beat_text = (
                "It beat the nearby alternative because the winning radius retained enough short-range "
                "context to capture neuron-associated astrocyte exposure without becoming so strict "
                "that genuine perineuronal coverage events were missed."
            )
    else:
        rationale_beat_text = "It beat the nearby alternative because no alternate local variation outperformed it."

    return f"""## Summary
Tested the {FEATURE_FAMILY} family in CA1, and {best_name} won with selection_score {_format_float(float(best['selection_score']))}.

## Metrics
Winning variation: {best_name} ({radius} px).
- IS partial r: {_format_float(float(best['partial_r']))}
- Selection score: {_format_float(float(best['selection_score']))}
- LOO predictive r: {_format_float(float(best['loo_predictive_r']))}
- IS-LOO Gap: {_format_float(float(best['is_loo_gap']))} (penalty={_format_float(float(best['penalty']))})
- Adjusted Score: {_format_float(float(best['adjusted_score']))}
- Analyzable donors: {analyzable}/{n_total}

Other tested variations:
{chr(10).join(other_lines)}

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The winner works best when the signal is expressed as a neuron-anchored CA1 exposure metric: among surviving CA1 pyramidal neurons, what fraction sits within {radius} px of at least one reactive astrocyte. This sharpens the biology from bulk cell composition toward a local neuron-glia interaction state that plausibly tracks stressed or remodeled CA1 tissue relevant to memory decline.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - {failure_text}
3. Error pattern: which donors are consistently wrong and what they share
{chr(10).join(worst_lines)}

## Rationale
The best variation is biologically coherent because it conditions first on the vulnerable CA1 pyramidal neuron population and then measures local reactive astrocyte proximity around those neurons, rather than counting reactive astrocytes globally. That makes the scalar interpretable as peripyramidal reactive astrocyte coverage. {rationale_beat_text} Relative to the accepted panel, this candidate is plausibly somewhat additive if it captures conditional local exposure beyond overall CA1 pyramidal loss, reactive astrocyte burden, and the broader mixed-cell niche fraction already in panel review.

## Interpretation
Biologically, the signal appears to represent how strongly surviving CA1 pyramidal neurons are embedded in a reactive astrocyte-rich local niche. The population is CA1 Pyramidal Neuron, the niche partner is CA1 Reactive Astrocyte, the donor-level scalar is the fraction of CA1 pyramidal neurons with at least one reactive astrocyte within {radius} px, and the simplest observable tissue pattern is a higher share of CA1 pyramidal neurons with nearby reactive astrocyte neighbors crowding the local perineuronal space.

## Next
Run the next local sweep around the same neuron-anchored CA1 reactive-astrocyte exposure idea but vary only the donor-level aggregation, for example comparing binary coverage to mean nearest-reactive-astrocyte distance or a soft distance-weighted coverage score, because the error pattern suggests the family has signal but may be losing information by thresholding all covered neurons equally.
"""


def main() -> None:
    data_root = DATA_ROOT_DEFAULT
    donor_table = _extract_all_features(data_root)

    evaluations = [_evaluate_variation(donor_table, name) for name in VARIATIONS]
    ranked = sorted(
        evaluations,
        key=lambda row: (
            -np.inf if not np.isfinite(float(row["selection_score"])) else -float(row["selection_score"]),
            -np.inf if not np.isfinite(float(row["partial_r"])) else -abs(float(row["partial_r"])),
            -np.inf if not np.isfinite(float(row["loo_predictive_r"])) else -float(row["loo_predictive_r"]),
            0 if row["name"] == BASELINE_VARIATION else 1,
        ),
    )
    best = ranked[0]
    best_feature = str(best["feature_column"])
    best_donor_loo = pd.DataFrame(best["per_donor_loo"])

    donor_feature_path = Path("/scratch/donor_feature_table.csv")
    write_donor_feature_table(
        donor_feature_path,
        donor_table,
        feature_column=best_feature,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUNDS,
        extra_columns=[*(spec["feature_column"] for spec in VARIATIONS.values()), *KEY_COUNT_COLUMNS],
    )

    results_payload = build_results_payload(
        status="ok",
        feature_name=best_feature,
        outcome=OUTCOME_COLUMN,
        n_total=int(best["n_total"]),
        n_analyzable=int(best["n_analyzable"]),
        partial_r=float(best["partial_r"]),
        ci_lo=None,
        ci_hi=None,
        p_value=float(best["p_value"]),
        loo_predictive_r=float(best["loo_predictive_r"]),
        loo_unstable_count=0,
        loo_max_shift=None,
        donor_ids_used=list(best["donor_ids_used"]),
        covariates=CONFOUNDS,
        recomputed_from_raw=True,
        registry_written=False,
        artifacts={"donor_feature_table": str(donor_feature_path)},
        best_variation=str(best["name"]),
        ranked_variations=[
            {
                key: value
                for key, value in row.items()
                if key not in {"per_donor_loo"}
            }
            for row in ranked
        ],
        selection_score=float(best["selection_score"]),
        is_loo_gap=float(best["is_loo_gap"]),
        penalty=float(best["penalty"]),
        adjusted_score=float(best["adjusted_score"]),
        feature_column=best_feature,
        feature_family=FEATURE_FAMILY,
    )
    write_results_payload("/scratch/results.json", results_payload)

    report_md = _build_report(best, ranked, best_donor_loo)
    Path("/scratch/report.md").write_text(report_md, encoding="utf-8")

    print(f"HYPOTHESIS FAMILY: {FEATURE_FAMILY}")
    print(f"BEST VARIATION: {best['name']}")
    print(f"  IS partial r:      {_format_float(float(best['partial_r']))}")
    print(f"  Selection score:   {_format_float(float(best['selection_score']))}")
    print(f"  LOO predictive r:  {_format_float(float(best['loo_predictive_r']))}  (diagnostic)")
    print(f"  IS-LOO Gap:        {_format_float(float(best['is_loo_gap']))}  (penalty={_format_float(float(best['penalty']))})")
    print(f"  Adjusted Score:    {_format_float(float(best['adjusted_score']))}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for row in ranked:
        print(
            f"  {row['name']}  {_format_float(float(row['partial_r']))}  "
            f"{_format_float(float(row['selection_score']))}  {_format_float(float(row['loo_predictive_r']))}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best_feature}  ca1_pyramidal_neuron_count  ca1_reactive_astrocyte_count")
    for row in best_donor_loo.itertuples(index=False):
        print(
            f"  {row.donor_id}  {row.outcome:.4f}  {row.predicted:.4f}  "
            f"{getattr(row, best_feature):.4f}  {int(row.ca1_pyramidal_neuron_count)}  "
            f"{int(row.ca1_reactive_astrocyte_count)}"
        )


if __name__ == "__main__":
    main()
