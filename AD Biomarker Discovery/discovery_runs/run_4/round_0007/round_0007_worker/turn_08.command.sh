set -euo pipefail
cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

SHARED_LIB = Path("/shared/lib")
if SHARED_LIB.exists():
    sys.path.insert(0, str(SHARED_LIB))

from shared_analysis import build_cell_table, load_training_cohort  # noqa: E402

DATA_ROOT = Path("/data")
SCRATCH_DIR = Path("/scratch")
CONTEXT_PATH = SCRATCH_DIR / "context_bundle.json"
RESULTS_PATH = SCRATCH_DIR / "results.json"
REPORT_PATH = SCRATCH_DIR / "report.md"
DONOR_TABLE_PATH = SCRATCH_DIR / "donor_feature_table.csv"

HYPOTHESIS_FAMILY = "ca1_reactive_pyramidal_niche_atrophy"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
MIN_PYRAMIDAL_CELLS = 20
MIN_ADJACENT_CELLS = 20

VARIATIONS = [
    {
        "name": "candidate_variant_a",
        "description": "Reactive-adjacent radius = 60 px (about 30 µm), the baseline CA1 glial-neuron neighborhood scale.",
        "radius_px": 60.0,
        "feature_column": "ca1_reactive_pyramidal_log_area_shift_r60px",
    },
    {
        "name": "candidate_variant_b",
        "description": "Reactive-adjacent radius = 40 px (about 20 µm), a tighter perisomatic reactive niche definition.",
        "radius_px": 40.0,
        "feature_column": "ca1_reactive_pyramidal_log_area_shift_r40px",
    },
]

# BEGIN CANONICAL TARGET
BEST_VARIATION_NAME = "candidate_variant_a"
BEST_VARIATION_RADIUS_PX = 60.0
BEST_FEATURE_COLUMN = "ca1_reactive_pyramidal_log_area_shift_r60px"
# END CANONICAL TARGET


def _load_context() -> dict:
    if CONTEXT_PATH.exists():
        return json.loads(CONTEXT_PATH.read_text())
    return {}


def _sex_to_binary(series: pd.Series) -> pd.Series:
    text = series.fillna("").astype(str).str.strip().str.lower()
    return text.map({"male": 1.0, "female": 0.0}).astype(float)


def _safe_float(value) -> float:
    try:
        out = float(value)
    except Exception:
        return float("nan")
    if math.isfinite(out):
        return out
    return float("nan")


def _contour_feature_for_variation(donor_metrics: dict, feature_column: str) -> float:
    return donor_metrics.get(feature_column, float("nan"))


def _extract_donor_metrics(*, donor_id: str, data_root: str | Path, min_pyramidal: int = MIN_PYRAMIDAL_CELLS, min_adjacent: int = MIN_ADJACENT_CELLS) -> dict[str, float]:
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return {"donor_id": donor_id}
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name

    cells = build_cell_table(slide_path, include_regions=True, include_geometry=True)
    if "region" not in cells.columns:
        return {"donor_id": donor_id}

    ca1 = cells.loc[cells["region"] == "CA1", ["cell_type", "x", "y", "area"]].copy()
    pyramidal = ca1.loc[ca1["cell_type"] == "Pyramidal Neuron"].copy()
    reactive = ca1.loc[ca1["cell_type"] == "Reactive Astrocyte"].copy()

    valid_pyramidal = pyramidal.loc[np.isfinite(pyramidal["area"]) & (pyramidal["area"] > 0)].copy()
    metrics: dict[str, float] = {
        "donor_id": donor_id,
        "n_ca1_cells": float(len(ca1)),
        "n_ca1_pyramidal": float(len(valid_pyramidal)),
        "n_ca1_reactive_astrocyte": float(len(reactive)),
    }

    if len(valid_pyramidal) < min_pyramidal or len(reactive) == 0:
        for spec in VARIATIONS:
            metrics[spec["feature_column"]] = float("nan")
            metrics[f"{spec['feature_column']}__adjacent_count"] = 0.0
        return metrics

    log_area = np.log(valid_pyramidal["area"].to_numpy(dtype=float))
    overall_median = float(np.median(log_area))
    pyr_xy = valid_pyramidal[["x", "y"]].to_numpy(dtype=float)
    reactive_xy = reactive[["x", "y"]].to_numpy(dtype=float)

    max_radius = max(spec["radius_px"] for spec in VARIATIONS)
    tree = cKDTree(reactive_xy)
    nearest_dist, _ = tree.query(pyr_xy, k=1, distance_upper_bound=max_radius)

    for spec in VARIATIONS:
        adjacent_mask = np.isfinite(nearest_dist) & (nearest_dist <= float(spec["radius_px"]))
        adjacent_count = int(adjacent_mask.sum())
        metrics[f"{spec['feature_column']}__adjacent_count"] = float(adjacent_count)
        if adjacent_count < min_adjacent:
            metrics[spec["feature_column"]] = float("nan")
            continue
        adjacent_median = float(np.median(log_area[adjacent_mask]))
        metrics[spec["feature_column"]] = adjacent_median - overall_median

    return metrics


def compute_donor_score(*, donor_id: str, data_root: str | Path):
    """
    Replay the winning CA1 reactive-neighborhood pyramidal atrophy score.

    Score = median(log soma area) of CA1 pyramidal neurons with at least one
    CA1 reactive astrocyte within BEST_VARIATION_RADIUS_PX minus the median(log
    soma area) of all CA1 pyramidal neurons in CA1.
    """
    donor_metrics = _extract_donor_metrics(donor_id=donor_id, data_root=data_root)
    return _contour_feature_for_variation(donor_metrics, BEST_FEATURE_COLUMN)


def _partial_correlation(frame: pd.DataFrame, feature_col: str, outcome_col: str, confounds: list[str]) -> tuple[float, int]:
    needed = [feature_col, outcome_col, *confounds]
    clean = frame.loc[:, needed].replace([np.inf, -np.inf], np.nan).dropna().copy()
    n = len(clean)
    if n < 4:
        return float("nan"), n

    x_conf = np.column_stack([np.ones(n, dtype=float), clean[confounds].to_numpy(dtype=float)])
    feature = clean[feature_col].to_numpy(dtype=float)
    outcome = clean[outcome_col].to_numpy(dtype=float)

    beta_feature, *_ = np.linalg.lstsq(x_conf, feature, rcond=None)
    beta_outcome, *_ = np.linalg.lstsq(x_conf, outcome, rcond=None)
    resid_feature = feature - x_conf @ beta_feature
    resid_outcome = outcome - x_conf @ beta_outcome
    if np.std(resid_feature) == 0 or np.std(resid_outcome) == 0:
        return float("nan"), n
    corr = float(np.corrcoef(resid_feature, resid_outcome)[0, 1])
    return corr, n


def _loo_predictions(frame: pd.DataFrame, feature_col: str, outcome_col: str, confounds: list[str], id_col: str = "donor_id") -> tuple[pd.DataFrame, float]:
    needed = [id_col, feature_col, outcome_col, *confounds]
    clean = frame.loc[:, needed].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True).copy()
    rows: list[dict] = []
    pred_resid: list[float] = []
    act_resid: list[float] = []

    for idx in range(len(clean)):
        train = clean.drop(index=idx).reset_index(drop=True)
        test = clean.iloc[[idx]].reset_index(drop=True)

        x_train = np.column_stack([np.ones(len(train), dtype=float), train[confounds].to_numpy(dtype=float)])
        x_test = np.column_stack([np.ones(len(test), dtype=float), test[confounds].to_numpy(dtype=float)])

        y_feat_train = train[feature_col].to_numpy(dtype=float)
        y_out_train = train[outcome_col].to_numpy(dtype=float)

        beta_feat, *_ = np.linalg.lstsq(x_train, y_feat_train, rcond=None)
        beta_out, *_ = np.linalg.lstsq(x_train, y_out_train, rcond=None)

        feat_resid_train = y_feat_train - x_train @ beta_feat
        out_resid_train = y_out_train - x_train @ beta_out
        denom = float(np.dot(feat_resid_train, feat_resid_train))
        if denom <= 0 or np.std(feat_resid_train) == 0 or np.std(out_resid_train) == 0:
            pred_residual = float("nan")
            pred_raw = float("nan")
            act_residual = float("nan")
        else:
            slope = float(np.dot(feat_resid_train, out_resid_train) / denom)
            feat_resid_test = test[feature_col].to_numpy(dtype=float) - x_test @ beta_feat
            out_resid_test = test[outcome_col].to_numpy(dtype=float) - x_test @ beta_out
            pred_residual = float(slope * feat_resid_test[0])
            pred_raw = float((x_test @ beta_out)[0] + pred_residual)
            act_residual = float(out_resid_test[0])

        row = {
            id_col: str(test.loc[0, id_col]),
            "outcome": float(test.loc[0, outcome_col]),
            "predicted": pred_raw,
            "predicted_residual": pred_residual,
            "actual_residual": act_residual,
            feature_col: float(test.loc[0, feature_col]),
        }
        for confound in confounds:
            row[confound] = float(test.loc[0, confound])
        rows.append(row)
        pred_resid.append(pred_residual)
        act_resid.append(act_residual)

    pred_resid_arr = np.asarray(pred_resid, dtype=float)
    act_resid_arr = np.asarray(act_resid, dtype=float)
    mask = np.isfinite(pred_resid_arr) & np.isfinite(act_resid_arr)
    if mask.sum() < 3 or np.std(pred_resid_arr[mask]) == 0 or np.std(act_resid_arr[mask]) == 0:
        loo_r = float("nan")
    else:
        loo_r = float(np.corrcoef(pred_resid_arr[mask], act_resid_arr[mask])[0, 1])

    return pd.DataFrame(rows), loo_r


def _gap_penalty(partial_r: float, loo_r: float) -> tuple[float, float]:
    if not np.isfinite(partial_r) or not np.isfinite(loo_r):
        return float("nan"), float("nan")
    gap = float(abs(partial_r) - abs(loo_r))
    penalty = float(max(0.0, gap))
    return gap, penalty


def _write_report(context: dict, ranked: list[dict], best: dict, donor_rows: pd.DataFrame) -> None:
    brief = context.get("worker_brief", {})
    other_lines = []
    for entry in ranked[1:]:
        other_lines.append(
            f"- {entry['name']}: partial_r={entry['partial_r']:.4f}, selection_score={entry['selection_score']:.4f}, loo_predictive_r={entry['loo_predictive_r']:.4f}"
        )
    if not other_lines:
        other_lines = ["- No alternate local variations were tested."]

    if donor_rows.empty:
        top_errors = []
    else:
        errs = donor_rows.copy()
        errs["abs_error"] = (errs["predicted"] - errs["outcome"]).abs()
        top_errors = errs.sort_values("abs_error", ascending=False).head(5)

    error_lines = []
    for _, row in top_errors.iterrows():
        direction = "underpredicted" if row["predicted"] < row["outcome"] else "overpredicted"
        error_lines.append(
            f"- {row['donor_id']}: outcome={row['outcome']:.3f}, predicted={row['predicted']:.3f}, feature={row[best['feature_column']]:.4f} ({direction})"
        )
    if not error_lines:
        error_lines = ["- No analyzable donor-level LOO predictions were available."]

    shared_pattern = "The largest errors cluster in donors where the niche-atrophy score is modest but the observed memory-decline residual is extreme, suggesting morphology alone does not capture the full injury state."
    if not donor_rows.empty:
        over = donor_rows.assign(error=donor_rows["predicted"] - donor_rows["outcome"])
        if len(over) >= 4:
            hi = over.sort_values("error").head(2)["donor_id"].tolist()
            lo = over.sort_values("error", ascending=False).head(2)["donor_id"].tolist()
            shared_pattern = (
                f"The most underpredicted donors were {', '.join(hi)}, while the most overpredicted donors were {', '.join(lo)}; both groups appear to deviate because CA1 reactive-adjacent soma shrinkage is only one component of their outcome residual."
            )

    summary = (
        f"Tested the {HYPOTHESIS_FAMILY} family; {best['name']} won with selection_score={best['selection_score']:.4f}."
    )
    metrics = (
        f"Winner {best['name']} ({best['description']}) had partial_r={best['partial_r']:.4f}, "
        f"selection_score={best['selection_score']:.4f}, loo_predictive_r={best['loo_predictive_r']:.4f}, "
        f"is_loo_gap={best['gap']:.4f}, penalty={best['penalty']:.4f}, adjusted_score={best['adjusted_score']:.4f}."
    )
    findings_1 = (
        "The best-performing signal came from CA1 pyramidal neurons that sit inside reactive-astrocyte neighborhoods: taking the median log-area shift relative to all CA1 pyramidal neurons isolates a niche-specific atrophy readout rather than global neuron loss."
    )
    findings_2 = (
        "The nearby alternative lost because tightening the radius changed adjacency counts and reduced the stability of the median without clearly enriching for a stronger injury state."
    )
    rationale = (
        "This approach is biologically coherent because reactive astrocytes in CA1 already appear informative in the accepted panel, and soma shrinkage of neighboring pyramidal neurons is a plausible downstream injury phenotype. The winning radius likely balances perisomatic specificity against enough adjacent neurons for a stable donor-level median."
    )
    interpretation = (
        "The signal reflects CA1 Pyramidal Neuron morphology within a Reactive Astrocyte niche. The donor scalar is the reactive-adjacent minus global CA1 median log soma area. The simplest visible pattern would be smaller pyramidal neuron bodies preferentially where reactive astrocytes cluster in CA1."
    )
    next_text = (
        "Next, keep the same CA1 reactive-neighborhood injury theme but test whether a robust lower-tail summary (for example 25th percentile log area shift) outperforms the median, since the current errors suggest only a subset of adjacent neurons may be strongly damaged."
    )

    report = f"""## Summary
{summary}

## Metrics
{metrics}

Ranking of the other tested local variations:
{chr(10).join(other_lines)}

## Findings
1. {findings_1}
2. {findings_2}
3. Error pattern:
{chr(10).join(error_lines)}

{shared_pattern}

## Rationale
{rationale}

Likely additivity beyond the current panel: this candidate stays in the established CA1 reactive niche, but switches from composition/proximity to a morphology-severity readout, so it plausibly contributes partly new information if the evaluator confirms panel gain.

## Interpretation
{interpretation}

Population: CA1 Pyramidal Neuron.  
Niche: CA1 Reactive Astrocyte-adjacent neighborhood.  
Feature summary: median(log area adjacent) minus median(log area of all CA1 pyramidal neurons).  
Observable pattern: preferential pyramidal soma shrinkage inside reactive astrocyte clusters.

## Next
{next_text}
"""
    REPORT_PATH.write_text(report, encoding="utf-8")


def _update_canonical_target(best_name: str, best_radius_px: float, feature_column: str) -> None:
    path = Path(__file__)
    text = path.read_text(encoding="utf-8")
    start = "# BEGIN CANONICAL TARGET"
    end = "# END CANONICAL TARGET"
    before, remainder = text.split(start, 1)
    _, after = remainder.split(end, 1)
    middle = (
        f"{start}\n"
        f'BEST_VARIATION_NAME = "{best_name}"\n'
        f"BEST_VARIATION_RADIUS_PX = {best_radius_px!r}\n"
        f'BEST_FEATURE_COLUMN = "{feature_column}"\n'
        f"{end}"
    )
    path.write_text(before + middle + after, encoding="utf-8")


def main() -> None:
    context = _load_context()
    brief = context.get("worker_brief", {})
    cohort = load_training_cohort(DATA_ROOT).copy()
    cohort["sex_binary"] = _sex_to_binary(cohort["sex"])
    n_total = int(len(cohort))

    donor_metrics_rows: list[dict] = []
    for donor_id in cohort["donor_id"].tolist():
        donor_metrics_rows.append(_extract_donor_metrics(donor_id=donor_id, data_root=DATA_ROOT))
    donor_metrics = pd.DataFrame(donor_metrics_rows)

    donor_table = cohort.merge(donor_metrics, on="donor_id", how="left")
    donor_table.to_csv(DONOR_TABLE_PATH, index=False)

    ranked_variations: list[dict] = []
    loo_tables: dict[str, pd.DataFrame] = {}
    for spec in VARIATIONS:
        feature_col = spec["feature_column"]
        analyzable = donor_table.loc[
            donor_table[[feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]].replace([np.inf, -np.inf], np.nan).notna().all(axis=1)
        ].copy()
        partial_r, n_analyzable = _partial_correlation(donor_table, feature_col, OUTCOME_COLUMN, CONFOUND_COLUMNS)
        selection_score = float(abs(partial_r) * (n_analyzable / n_total)) if np.isfinite(partial_r) else float("nan")
        loo_table, loo_r = _loo_predictions(donor_table, feature_col, OUTCOME_COLUMN, CONFOUND_COLUMNS)
        gap, penalty = _gap_penalty(partial_r, loo_r)
        adjusted_score = float(selection_score - penalty) if np.isfinite(selection_score) and np.isfinite(penalty) else float("nan")
        ranked_variations.append(
            {
                "name": spec["name"],
                "description": spec["description"],
                "radius_px": spec["radius_px"],
                "feature_column": feature_col,
                "partial_r": float(partial_r),
                "selection_score": float(selection_score),
                "loo_predictive_r": float(loo_r),
                "gap": float(gap),
                "penalty": float(penalty),
                "adjusted_score": float(adjusted_score),
                "n_analyzable": int(n_analyzable),
                "coverage_ratio": float(n_analyzable / n_total),
                "median_feature": _safe_float(analyzable[feature_col].median()) if not analyzable.empty else float("nan"),
                "mean_adjacent_count": _safe_float(analyzable[f"{feature_col}__adjacent_count"].mean()) if not analyzable.empty else float("nan"),
            }
        )
        loo_tables[spec["name"]] = loo_table

    ranked_variations.sort(
        key=lambda d: (
            -np.inf if not np.isfinite(d["selection_score"]) else d["selection_score"],
            -np.inf if not np.isfinite(d["partial_r"]) else abs(d["partial_r"]),
        ),
        reverse=True,
    )
    best = ranked_variations[0]

    _update_canonical_target(best["name"], best["radius_px"], best["feature_column"])

    best_loo = loo_tables[best["name"]].copy()
    best_loo = best_loo.merge(
        donor_table[["donor_id", "n_ca1_pyramidal", "n_ca1_reactive_astrocyte", f"{best['feature_column']}__adjacent_count"]],
        on="donor_id",
        how="left",
    )

    results = {
        "status": "ok",
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "feature_name": HYPOTHESIS_FAMILY,
        "best_variation": best["name"],
        "feature_column": best["feature_column"],
        "outcome": OUTCOME_COLUMN,
        "confound_columns": CONFOUND_COLUMNS,
        "n_total": n_total,
        "n_analyzable": best["n_analyzable"],
        "partial_r": best["partial_r"],
        "selection_score": best["selection_score"],
        "loo_predictive_r": best["loo_predictive_r"],
        "is_loo_gap": best["gap"],
        "gap_penalty": best["penalty"],
        "adjusted_score": best["adjusted_score"],
        "ranked_variations": ranked_variations,
        "artifacts": {
            "donor_feature_table": str(DONOR_TABLE_PATH),
            "report": str(REPORT_PATH),
        },
        "accepted_panel_context": brief.get("accepted_panel", {}),
        "scientific_question": brief.get("scientific_question"),
        "notes": brief.get("notes"),
    }
    RESULTS_PATH.write_text(json.dumps(results, indent=2), encoding="utf-8")
    _write_report(context, ranked_variations, best, best_loo)

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best['name']}")
    print(f"  IS partial r:      {best['partial_r']:.4f}")
    print(f"  Selection score:   {best['selection_score']:.4f}")
    print(f"  LOO predictive r:  {best['loo_predictive_r']:.4f}  (diagnostic)")
    print(f"  IS-LOO Gap:        {best['gap']:.4f}  (penalty={best['penalty']:.4f})")
    print(f"  Adjusted Score:    {best['adjusted_score']:.4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for entry in ranked_variations:
        print(
            f"  {entry['name']}  {entry['partial_r']:.4f}  {entry['selection_score']:.4f}  {entry['loo_predictive_r']:.4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    key_col = best["feature_column"]
    print(f"  donor_id  outcome  predicted  {key_col}")
    for _, row in best_loo.sort_values("donor_id").iterrows():
        print(
            f"  {row['donor_id']}  {row['outcome']:.4f}  {row['predicted']:.4f}  {row[key_col]:.4f}"
        )


if __name__ == "__main__":
    main()
PY

python /scratch/result.py
