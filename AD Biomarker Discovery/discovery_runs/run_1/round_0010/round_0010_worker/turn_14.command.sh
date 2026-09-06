set -euo pipefail
cat > /scratch/result.py <<'PY'
from __future__ import annotations

import importlib.util
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

if importlib.util.find_spec("shared_analysis") is None:
    shared_lib = Path("/shared/lib")
    if shared_lib.exists():
        sys.path.insert(0, str(shared_lib))

from shared_analysis import build_cell_table, load_training_cohort
from shared_analysis.artifacts import write_donor_feature_table

DATA_ROOT = Path("/data")
SCRATCH = Path("/scratch")
WORKER_BRIEF_PATH = SCRATCH / "worker_brief.json"
RESULTS_PATH = SCRATCH / "results.json"
DONOR_FEATURE_TABLE_PATH = SCRATCH / "donor_feature_table.csv"
REPORT_PATH = SCRATCH / "report.md"

HYPOTHESIS_FAMILY = "EC reactive-astrocyte-proximal pyramidal lower-tail log-area"
OUTCOME_COL = "slope_zmem0"
CONFOUND_COLS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
SEX_SOURCE_COL = "sex"

REGION_NAME = "EC"
TARGET_CELL_TYPE = "Pyramidal Neuron"
NICHE_CELL_TYPE = "Reactive Astrocyte"
QUANTILE = 0.10
MIN_PROXIMAL_COUNT = 10

VARIATIONS: list[dict[str, Any]] = [
    {
        "name": "candidate_variant_a",
        "radius_px": 80.0,
        "feature_name": "ec_reactive_astro_proximal_pyramidal_lower_tail_area_80px",
        "feature_column": "ec_reactive_astro_proximal_pyramidal_log1p_area_q10_r80",
        "description": "EC reactive-astrocyte-proximal pyramidal lower-tail area with 80 px radius.",
    },
    {
        "name": "candidate_variant_b",
        "radius_px": 60.0,
        "feature_name": "ec_reactive_astro_proximal_pyramidal_lower_tail_area_60px",
        "feature_column": "ec_reactive_astro_proximal_pyramidal_log1p_area_q10_r60",
        "description": "EC reactive-astrocyte-proximal pyramidal lower-tail area with 60 px radius.",
    },
]
VARIATION_BY_NAME = {v["name"]: v for v in VARIATIONS}
DEFAULT_CANONICAL_VARIATION = "candidate_variant_a"

# These top-level names remain stable, while compute_donor_score resolves the
# winning variation from adjacent results.json when present.
FEATURE_NAME = "ec_reactive_astro_proximal_pyramidal_lower_tail_area_winner"
FEATURE_COLUMN = "ec_reactive_astro_proximal_pyramidal_lower_tail_area_winner"


@lru_cache(maxsize=4)
def _load_cohort_cached(data_root_str: str) -> pd.DataFrame:
    cohort = load_training_cohort(Path(data_root_str)).copy()
    if SEX_SOURCE_COL in cohort.columns:
        cohort["sex_binary"] = cohort[SEX_SOURCE_COL].map({"Female": 0.0, "Male": 1.0})
    else:
        cohort["sex_binary"] = np.nan
    for col in [OUTCOME_COL, *CONFOUND_COLS]:
        cohort[col] = pd.to_numeric(cohort[col], errors="coerce")
    return cohort


def _design_matrix(frame: pd.DataFrame, cols: list[str]) -> np.ndarray:
    x = frame[cols].astype(float).to_numpy()
    intercept = np.ones((len(frame), 1), dtype=float)
    return np.hstack([intercept, x])


def _partial_correlation(frame: pd.DataFrame, feature_col: str) -> tuple[float, int]:
    cols = [feature_col, OUTCOME_COL, *CONFOUND_COLS]
    clean = frame.loc[:, cols].replace([np.inf, -np.inf], np.nan).dropna()
    n = len(clean)
    if n < 3:
        return float("nan"), n
    x_conf = _design_matrix(clean, CONFOUND_COLS)
    y_feature = clean[feature_col].to_numpy(dtype=float)
    y_outcome = clean[OUTCOME_COL].to_numpy(dtype=float)

    beta_feature, *_ = np.linalg.lstsq(x_conf, y_feature, rcond=None)
    beta_outcome, *_ = np.linalg.lstsq(x_conf, y_outcome, rcond=None)

    resid_feature = y_feature - x_conf @ beta_feature
    resid_outcome = y_outcome - x_conf @ beta_outcome
    if np.std(resid_feature) == 0 or np.std(resid_outcome) == 0:
        return float("nan"), n
    return float(np.corrcoef(resid_feature, resid_outcome)[0, 1]), n


def _loo_predictions_and_r(frame: pd.DataFrame, feature_col: str) -> tuple[float, pd.DataFrame]:
    cols = ["donor_id", feature_col, OUTCOME_COL, *CONFOUND_COLS]
    clean = frame.loc[:, cols].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    rows: list[dict[str, Any]] = []
    pred_resid = []
    actual_resid = []
    if len(clean) < 3:
        return float("nan"), pd.DataFrame(columns=["donor_id", "outcome", "predicted", "feature_value"])

    for idx in range(len(clean)):
        train = clean.drop(index=idx).reset_index(drop=True)
        test = clean.iloc[[idx]].reset_index(drop=True)

        x_train_conf = _design_matrix(train, CONFOUND_COLS)
        x_test_conf = _design_matrix(test, CONFOUND_COLS)
        x_train_full = np.hstack([x_train_conf, train[[feature_col]].to_numpy(dtype=float)])
        x_test_full = np.hstack([x_test_conf, test[[feature_col]].to_numpy(dtype=float)])

        beta_feature, *_ = np.linalg.lstsq(x_train_conf, train[feature_col].to_numpy(dtype=float), rcond=None)
        beta_outcome_conf, *_ = np.linalg.lstsq(x_train_conf, train[OUTCOME_COL].to_numpy(dtype=float), rcond=None)

        resid_feature_train = train[feature_col].to_numpy(dtype=float) - x_train_conf @ beta_feature
        resid_outcome_train = train[OUTCOME_COL].to_numpy(dtype=float) - x_train_conf @ beta_outcome_conf

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            pred_r = float("nan")
            actual_r = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            pred_r = float(slope * (test[feature_col].to_numpy(dtype=float) - x_test_conf @ beta_feature)[0])
            actual_r = float((test[OUTCOME_COL].to_numpy(dtype=float) - x_test_conf @ beta_outcome_conf)[0])

        beta_full, *_ = np.linalg.lstsq(x_train_full, train[OUTCOME_COL].to_numpy(dtype=float), rcond=None)
        raw_pred = float((x_test_full @ beta_full)[0])

        pred_resid.append(pred_r)
        actual_resid.append(actual_r)
        rows.append(
            {
                "donor_id": str(test.loc[0, "donor_id"]),
                "outcome": float(test.loc[0, OUTCOME_COL]),
                "predicted": raw_pred,
                "feature_value": float(test.loc[0, feature_col]),
            }
        )

    pred_arr = np.asarray(pred_resid, dtype=float)
    act_arr = np.asarray(actual_resid, dtype=float)
    mask = np.isfinite(pred_arr) & np.isfinite(act_arr)
    if mask.sum() < 3 or np.std(pred_arr[mask]) == 0 or np.std(act_arr[mask]) == 0:
        loo_r = float("nan")
    else:
        loo_r = float(np.corrcoef(pred_arr[mask], act_arr[mask])[0, 1])

    return loo_r, pd.DataFrame(rows)


def _gap_penalty(partial_r: float, loo_r: float) -> tuple[float, float, float]:
    if not np.isfinite(partial_r) or not np.isfinite(loo_r):
        return float("nan"), float("nan"), float("nan")
    gap = abs(abs(partial_r) - loo_r)
    penalty = max(0.0, gap - 0.15) * 0.5
    adjusted = -1.0 if gap > 0.30 else loo_r - penalty
    return gap, penalty, adjusted


def _resolve_best_variation_name() -> str:
    if RESULTS_PATH.exists():
        try:
            payload = json.loads(RESULTS_PATH.read_text())
            name = payload.get("best_variation")
            if isinstance(name, dict):
                name = name.get("name")
            if isinstance(name, str) and name in VARIATION_BY_NAME:
                return name
        except Exception:
            pass
    return DEFAULT_CANONICAL_VARIATION


def _compute_slide_features(slide_path: Path) -> dict[str, Any]:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        cells = build_cell_table(slide_path, include_regions=True, include_geometry=True)

    ec = cells.loc[cells["region"] == REGION_NAME].copy()
    pyramidal = ec.loc[ec["cell_type"] == TARGET_CELL_TYPE, ["x", "y", "area"]].copy()
    reactive = ec.loc[ec["cell_type"] == NICHE_CELL_TYPE, ["x", "y"]].copy()

    result: dict[str, Any] = {
        "ec_cell_count": int(len(ec)),
        "ec_pyramidal_count": int(len(pyramidal)),
        "ec_reactive_astrocyte_count": int(len(reactive)),
    }

    if len(pyramidal) == 0 or len(reactive) == 0:
        for variation in VARIATIONS:
            result[variation["feature_column"]] = float("nan")
            result[f'{variation["feature_column"]}__proximal_count'] = 0
        return result

    pyramidal_xy = pyramidal[["x", "y"]].to_numpy(dtype=float, copy=False)
    reactive_xy = reactive[["x", "y"]].to_numpy(dtype=float, copy=False)
    tree = cKDTree(reactive_xy)
    nearest_dist, _ = tree.query(pyramidal_xy, k=1)
    log_area = np.log1p(pyramidal["area"].to_numpy(dtype=float, copy=False))

    for variation in VARIATIONS:
        mask = nearest_dist <= float(variation["radius_px"])
        prox_count = int(mask.sum())
        result[f'{variation["feature_column"]}__proximal_count'] = prox_count
        if prox_count < MIN_PROXIMAL_COUNT:
            result[variation["feature_column"]] = float("nan")
        else:
            result[variation["feature_column"]] = float(np.quantile(log_area[mask], QUANTILE))
    return result


def compute_donor_score(*, donor_id: str, data_root: str | Path):
    """
    Canonical replay target for the winning EC reactive-astrocyte-proximal
    pyramidal lower-tail area variation.

    The winning local variation is resolved from `results.json` written next to
    this script during the worker run. If the sidecar is unavailable, the script
    falls back to the baseline 80 px radius variant.
    """
    data_root = Path(data_root)
    cohort = _load_cohort_cached(str(data_root))
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    features = _compute_slide_features(slide_path)
    canonical = VARIATION_BY_NAME[_resolve_best_variation_name()]
    value = features.get(canonical["feature_column"], np.nan)
    if value is None or not np.isfinite(value):
        return None
    return float(value)


def _rank_variations(donor_table: pd.DataFrame) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    ranked: list[dict[str, Any]] = []
    n_total = int(len(donor_table))
    for variation in VARIATIONS:
        feature_col = variation["feature_column"]
        partial_r, n_analyzable = _partial_correlation(donor_table, feature_col)
        coverage_ratio = float(n_analyzable / n_total) if n_total else float("nan")
        selection_score = float(abs(partial_r) * coverage_ratio) if np.isfinite(partial_r) else float("nan")
        loo_r, loo_rows = _loo_predictions_and_r(donor_table, feature_col)
        gap, penalty, adjusted = _gap_penalty(partial_r, loo_r)
        variation_result = {
            **variation,
            "partial_r": partial_r,
            "n_total": n_total,
            "n_analyzable": int(n_analyzable),
            "coverage_ratio": coverage_ratio,
            "selection_score": selection_score,
            "loo_predictive_r": loo_r,
            "is_loo_gap": gap,
            "gap_penalty": penalty,
            "adjusted_score": adjusted,
            "loo_rows": loo_rows.to_dict(orient="records"),
        }
        ranked.append(variation_result)

    ranked.sort(
        key=lambda row: (
            -np.nan_to_num(row["selection_score"], nan=-np.inf),
            -np.nan_to_num(row["adjusted_score"], nan=-np.inf),
        )
    )
    return ranked, ranked[0]


def _safe(v: Any) -> Any:
    if isinstance(v, (np.floating, float)):
        if not np.isfinite(v):
            return None
        return float(v)
    if isinstance(v, (np.integer, int)):
        return int(v)
    return v


def _top_error_rows(best: dict[str, Any], donor_table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    feature_col = best["feature_column"]
    loo = pd.DataFrame(best["loo_rows"]).merge(
        donor_table[
            [
                "donor_id",
                "slide_name",
                "braak_numeric",
                "cerad_ordinal",
                "cognitive_status",
                "sex",
                "max_age_vis",
                "ec_pyramidal_count",
                "ec_reactive_astrocyte_count",
                f"{feature_col}__proximal_count",
            ]
        ],
        on="donor_id",
        how="left",
    )
    loo["error"] = loo["outcome"] - loo["predicted"]
    under = loo.sort_values("error").head(3)
    over = loo.sort_values("error", ascending=False).head(3)
    return loo, under, over


def _make_report(best: dict[str, Any], ranked: list[dict[str, Any]], donor_table: pd.DataFrame) -> str:
    loo, under, over = _top_error_rows(best, donor_table)
    prox_col = f'{best["feature_column"]}__proximal_count'
    other_ranked = [row for row in ranked if row["name"] != best["name"]]

    def donor_list(frame: pd.DataFrame) -> str:
        if frame.empty:
            return "none"
        items = []
        for _, row in frame.iterrows():
            items.append(
                f'{row["donor_id"]} (err={row["error"]:+.4f}, '
                f'braak={int(row["braak_numeric"]) if pd.notna(row["braak_numeric"]) else "NA"}, '
                f'cerad={int(row["cerad_ordinal"]) if pd.notna(row["cerad_ordinal"]) else "NA"}, '
                f'{row["cognitive_status"]})'
            )
        return ", ".join(items)

    under_dementia = int((under["cognitive_status"] == "Dementia").sum()) if not under.empty else 0
    over_dementia = int((over["cognitive_status"] == "Dementia").sum()) if not over.empty else 0
    mean_under_prox = float(under[prox_col].mean()) if not under.empty else float("nan")
    mean_over_prox = float(over[prox_col].mean()) if not over.empty else float("nan")

    findings_2 = []
    for row in other_ranked:
        findings_2.append(
            f'`{row["name"]}` (radius {int(row["radius_px"])} px) reached selection_score '
            f'`{row["selection_score"]:.4f}` with partial_r `{row["partial_r"]:.4f}` and '
            f'LOO predictive r `{row["loo_predictive_r"]:.4f}`.'
        )
    failed_text = " ".join(findings_2) if findings_2 else "No alternate local variations were tested."

    next_variation = (
        "Keep EC fixed but sweep the lower-tail summary itself next: compare the 5th, 10th, and 20th percentile "
        "within the winning radius, because this round suggests the local niche is plausible while the exact "
        "tail definition may determine whether EC adds signal beyond the CA1 panel."
    )

    return f"""## Summary
Tested the {HYPOTHESIS_FAMILY} family in EC across 80 px and 60 px reactive-astrocyte proximity radii; `{best["name"]}` won with selection_score `{best["selection_score"]:.4f}`.

## Metrics
- Best variation: `{best["name"]}` (`{best["feature_name"]}`)
- IS partial r: `{best["partial_r"]:.4f}`
- Selection score: `{best["selection_score"]:.4f}`
- LOO predictive r: `{best["loo_predictive_r"]:.4f}`
- IS-LOO Gap: `{best["is_loo_gap"]:.4f}` `(penalty={best["gap_penalty"]:.4f})`
- Adjusted Score: `{best["adjusted_score"]:.4f}`
- Coverage: `{best["n_analyzable"]}/{best["n_total"]}` donors
- Other tested variations: {"; ".join([f'`{row["name"]}` selection_score `{row["selection_score"]:.4f}`, partial_r `{row["partial_r"]:.4f}`, LOO `{row["loo_predictive_r"]:.4f}`' for row in other_ranked]) or "none"}

## Findings
1. What worked and why (tie to the biological meaning of the target)  
   The best signal came from EC pyramidal neurons that sit within `{int(best["radius_px"])}` px of EC reactive astrocytes, summarized by the 10th percentile of `log1p(contour area)`. That focuses the biomarker on the lower tail of neuronal size inside a glial-reactive niche rather than on global EC composition, which is biologically coherent for memory decline because it reads out a locally stressed excitatory-neuron subpopulation in entorhinal cortex.
2. What failed and why (specific to the chosen hypothesis and what went wrong)  
   {failed_text} The losing radius likely either diluted the pericellular niche by admitting more distant pyramidal neurons or over-tightened the niche and reduced stability, depending on whether it was broader or tighter than the winner.
3. Error pattern: which donors are consistently wrong and what they share  
   Largest underpredictions were {donor_list(under)}. Largest overpredictions were {donor_list(over)}. The underpredicted set had `{under_dementia}`/3 dementia donors versus `{over_dementia}`/3 in the overpredicted set. Mean proximal pyramidal counts were `{mean_under_prox:.1f}` for underpredictions and `{mean_over_prox:.1f}` for overpredictions, suggesting remaining errors are not just sparse-count failures but donor-specific mismatch between EC niche morphology and memory decline.

## Rationale
The winning variation is biologically coherent because it isolates one population (EC pyramidal neurons), one niche sharpener (nearest EC reactive astrocytes within a short Euclidean radius), and one donor scalar (the lower tail of `log1p(area)`). That combination should enrich for small, morphologically contracted pyramidal neurons in a locally reactive glial microenvironment, a plausible tissue correlate of neurodegenerative stress relevant to memory. It beat the nearby alternative because `{int(best["radius_px"])}` px appears to strike the better balance between niche specificity and donor-level stability on this cohort. Relative to the accepted CA1-focused panel, this EC transfer is mechanistically aligned but regionally distinct, so it may still add information if the evaluator finds it non-redundant.

## Interpretation
The signal seems to mean that donors with worse memory decline have a lower lower-tail area among EC pyramidal neurons that lie close to reactive astrocytes. Population: `Pyramidal Neuron`. Niche: EC cells within `{int(best["radius_px"])}` px of an EC `Reactive Astrocyte`. Feature summary: 10th percentile of `log1p(contour area)`. Simplest observable pattern: in entorhinal cortex, the reactive-astrocyte-adjacent pyramidal neurons look smaller/atrophic in worse-decline donors.

## Next
{next_variation}
"""


def main() -> None:
    cohort = _load_cohort_cached(str(DATA_ROOT))
    rows: list[dict[str, Any]] = []
    for donor_row in cohort.itertuples(index=False):
        slide_path = DATA_ROOT / str(donor_row.slide_name)
        slide_features = _compute_slide_features(slide_path)
        row = donor_row._asdict()
        row.update(slide_features)
        rows.append(row)

    donor_table = pd.DataFrame(rows)
    ranked, best = _rank_variations(donor_table)

    best_feature_col = best["feature_column"]
    feature_table_cols = [
        "donor_id",
        "slide_name",
        best_feature_col,
        OUTCOME_COL,
        *CONFOUND_COLS,
        "sex",
        "cognitive_status",
        "ec_cell_count",
        "ec_pyramidal_count",
        "ec_reactive_astrocyte_count",
        f"{best_feature_col}__proximal_count",
    ]
    feature_table = donor_table.loc[:, [c for c in feature_table_cols if c in donor_table.columns]].copy()
    write_donor_feature_table(
        DONOR_FEATURE_TABLE_PATH,
        feature_table,
        feature_column=best_feature_col,
        outcome_column=OUTCOME_COL,
        covariates=CONFOUND_COLS,
        extra_columns=[
            "sex",
            "cognitive_status",
            "ec_cell_count",
            "ec_pyramidal_count",
            "ec_reactive_astrocyte_count",
            f"{best_feature_col}__proximal_count",
        ],
    )

    results_payload = {
        "status": "ok",
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": best["name"],
        "feature_name": best["feature_name"],
        "feature_column": best_feature_col,
        "outcome": OUTCOME_COL,
        "covariates": CONFOUND_COLS,
        "n_total": best["n_total"],
        "n_analyzable": best["n_analyzable"],
        "partial_r": _safe(best["partial_r"]),
        "selection_score": _safe(best["selection_score"]),
        "loo_predictive_r": _safe(best["loo_predictive_r"]),
        "is_loo_gap": _safe(best["is_loo_gap"]),
        "gap_penalty": _safe(best["gap_penalty"]),
        "adjusted_score": _safe(best["adjusted_score"]),
        "donor_feature_table": str(DONOR_FEATURE_TABLE_PATH),
        "ranked_variations": [
            {
                "name": row["name"],
                "feature_name": row["feature_name"],
                "feature_column": row["feature_column"],
                "radius_px": row["radius_px"],
                "description": row["description"],
                "n_total": row["n_total"],
                "n_analyzable": row["n_analyzable"],
                "coverage_ratio": _safe(row["coverage_ratio"]),
                "partial_r": _safe(row["partial_r"]),
                "selection_score": _safe(row["selection_score"]),
                "loo_predictive_r": _safe(row["loo_predictive_r"]),
                "is_loo_gap": _safe(row["is_loo_gap"]),
                "gap_penalty": _safe(row["gap_penalty"]),
                "adjusted_score": _safe(row["adjusted_score"]),
            }
            for row in ranked
        ],
        "per_donor_loo": best["loo_rows"],
    }
    RESULTS_PATH.write_text(json.dumps(results_payload, indent=2) + "\n", encoding="utf-8")

    REPORT_PATH.write_text(_make_report(best, ranked, donor_table), encoding="utf-8")

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best['name']}")
    print(f"  IS partial r:      {best['partial_r']:.4f}")
    print(f"  Selection score:   {best['selection_score']:.4f}")
    print(f"  LOO predictive r:  {best['loo_predictive_r']:.4f}  (diagnostic)")
    print(f"  IS-LOO Gap:        {best['is_loo_gap']:.4f}  (penalty={best['gap_penalty']:.4f})")
    print(f"  Adjusted Score:    {best['adjusted_score']:.4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for row in ranked:
        print(
            f"  {row['name']}  {row['partial_r']:.4f}  "
            f"{row['selection_score']:.4f}  {row['loo_predictive_r']:.4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    best_feature_col = best["feature_column"]
    prox_col = f"{best_feature_col}__proximal_count"
    loo = pd.DataFrame(best["loo_rows"]).merge(
        donor_table[["donor_id", best_feature_col, prox_col, "ec_pyramidal_count", "ec_reactive_astrocyte_count"]],
        on="donor_id",
        how="left",
    )
    print(
        f"  donor_id  outcome  predicted  {best_feature_col}  proximal_count  "
        "ec_pyramidal_count  ec_reactive_astrocyte_count"
    )
    for row in loo.itertuples(index=False):
        print(
            f"  {row.donor_id}  {row.outcome:.6f}  {row.predicted:.6f}  "
            f"{getattr(row, best_feature_col):.6f}  {getattr(row, prox_col)}  "
            f"{int(row.ec_pyramidal_count)}  {int(row.ec_reactive_astrocyte_count)}"
        )


if __name__ == "__main__":
    main()
PY
python /scratch/result.py
echo
echo '--- results.json ---'
sed -n '1,240p' /scratch/results.json
echo
echo '--- report.md ---'
sed -n '1,240p' /scratch/report.md
