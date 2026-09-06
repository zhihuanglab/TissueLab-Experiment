set -euo pipefail
cat > /scratch/result.py <<'PY'
from __future__ import annotations

import gc
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/shared/lib")

from shared_analysis.artifacts import build_results_payload, write_donor_feature_table, write_results_payload
from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort
from shared_analysis.stats import leave_one_out_summary, partial_correlation


FAMILY_NAME = "ca1_reactive_astrocyte_hypertrophy"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUNDS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
DEFAULT_VARIATION = "median_log_area_shift"
MIN_CELLS_PER_GROUP = 20

VARIATIONS = [
    {
        "name": "median_log_area_shift",
        "description": "median(log contour area) of CA1 Reactive Astrocytes minus median(log contour area) of CA1 Astrocytes",
        "quantile": 0.50,
        "feature_column": "ca1_reactive_astrocyte_hypertrophy_median_log_area_shift",
    },
    {
        "name": "upper_tail_log_area_shift",
        "description": "75th percentile(log contour area) of CA1 Reactive Astrocytes minus 75th percentile(log contour area) of CA1 Astrocytes",
        "quantile": 0.75,
        "feature_column": "ca1_reactive_astrocyte_hypertrophy_upper_tail_log_area_shift",
    },
]
VARIATION_BY_NAME = {row["name"]: row for row in VARIATIONS}

FEATURE_NAME = FAMILY_NAME
FEATURE_COLUMN = VARIATION_BY_NAME[DEFAULT_VARIATION]["feature_column"]

warnings.filterwarnings("ignore", message="Object at .* is not recognized as a component of a Zarr hierarchy.")
warnings.filterwarnings("ignore", message="invalid value encountered in divide", category=RuntimeWarning)


def _paths(base_dir: str | Path | None = None) -> dict[str, Path]:
    root = Path(base_dir) if base_dir is not None else Path(__file__).resolve().parent
    return {
        "root": root,
        "results": root / "results.json",
        "table": root / "donor_feature_table.csv",
        "report": root / "report.md",
        "state": root / "result_state.json",
    }


def _load_state(base_dir: str | Path | None = None) -> dict[str, object]:
    path = _paths(base_dir)["state"]
    if path.exists():
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}


def _canonical_variation(base_dir: str | Path | None = None) -> str:
    state = _load_state(base_dir)
    variation = state.get("best_variation")
    if isinstance(variation, str) and variation in VARIATION_BY_NAME:
        return variation
    return DEFAULT_VARIATION


def _feature_column_for(variation_name: str) -> str:
    return str(VARIATION_BY_NAME[variation_name]["feature_column"])


def _selection_score(partial_r: float | None, n_analyzable: int | None, n_total: int | None) -> float:
    if partial_r is None or not np.isfinite(partial_r):
        return float("nan")
    if n_analyzable is None or n_total is None or n_total <= 0:
        return float("nan")
    coverage = max(0.0, min(1.0, float(n_analyzable) / float(n_total)))
    return float(abs(float(partial_r)) * coverage)


def _adjusted_metrics(partial_r: float | None, loo_predictive_r: float | None) -> tuple[float, float, float]:
    is_r = abs(float(partial_r)) if partial_r is not None and np.isfinite(partial_r) else 0.0
    loo_r = abs(float(loo_predictive_r)) if loo_predictive_r is not None and np.isfinite(loo_predictive_r) else 0.0
    gap = is_r - loo_r
    if gap > 0.30:
        return float(gap), float(gap), -1.0
    penalty = max(0.0, gap - 0.15) * 0.5
    return float(gap), float(penalty), float(loo_r - penalty)


def _load_eval_cohort(data_root: str | Path) -> pd.DataFrame:
    data_root = Path(data_root)
    training = data_root / "training_cohort.csv"
    test = data_root / "test_cohort.csv"
    if training.exists():
        return load_training_cohort(data_root)
    if test.exists():
        return pd.read_csv(test)
    return load_training_cohort(data_root)


def _cohort_with_covariates(data_root: str | Path) -> pd.DataFrame:
    cohort = _load_eval_cohort(data_root).copy()
    if "sex_binary" not in cohort.columns:
        cohort["sex_binary"] = cohort["sex"].map({"Female": 1.0, "Male": 0.0}).astype(float)
    return cohort


def _safe_quantile(values: np.ndarray, q: float) -> float:
    if values.size == 0:
        return float("nan")
    return float(np.quantile(values, q))


def _extract_donor_summary(slide_path: str | Path) -> dict[str, float]:
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=True)
    mask = (
        (cells["region"] == "CA1")
        & (cells["cell_type"].isin(["Astrocyte", "Reactive Astrocyte"]))
        & np.isfinite(cells["area"].to_numpy(dtype=float))
        & (cells["area"].to_numpy(dtype=float) > 0)
    )
    subset = cells.loc[mask, ["cell_type", "area"]].copy()
    del cells
    gc.collect()

    result: dict[str, float] = {
        "ca1_astrocyte_count": 0.0,
        "ca1_reactive_astrocyte_count": 0.0,
        "ca1_total_astroglia_count": 0.0,
    }
    for variation in VARIATIONS:
        result[variation["feature_column"]] = float("nan")

    if subset.empty:
        return result

    subset["log_area"] = np.log(subset["area"].to_numpy(dtype=float))
    astro = subset.loc[subset["cell_type"] == "Astrocyte", "log_area"].to_numpy(dtype=float)
    reactive = subset.loc[subset["cell_type"] == "Reactive Astrocyte", "log_area"].to_numpy(dtype=float)

    result["ca1_astrocyte_count"] = float(len(astro))
    result["ca1_reactive_astrocyte_count"] = float(len(reactive))
    result["ca1_total_astroglia_count"] = float(len(astro) + len(reactive))

    if len(astro) < MIN_CELLS_PER_GROUP or len(reactive) < MIN_CELLS_PER_GROUP:
        return result

    for variation in VARIATIONS:
        q = float(variation["quantile"])
        result[variation["feature_column"]] = _safe_quantile(reactive, q) - _safe_quantile(astro, q)
    return result


def _design_matrix(frame: pd.DataFrame) -> np.ndarray:
    conf = frame.to_numpy(dtype=float)
    return np.column_stack([np.ones(len(frame), dtype=float), conf])


def _loo_predictions(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
    id_col: str = "donor_id",
) -> tuple[float, pd.DataFrame]:
    columns = [id_col, feature_col, outcome_col, *confounds]
    frame = df.loc[:, columns].dropna().reset_index(drop=True).copy()
    if frame.empty:
        return float("nan"), pd.DataFrame(columns=[id_col, outcome_col, "predicted", feature_col])

    rows: list[dict[str, float | str]] = []
    pred_resid: list[float] = []
    act_resid: list[float] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train_conf = _design_matrix(train[confounds])
        x_test_conf = _design_matrix(test[confounds])

        y_feature_train = train[feature_col].to_numpy(dtype=float)
        y_outcome_train = train[outcome_col].to_numpy(dtype=float)
        y_feature_test = test[feature_col].to_numpy(dtype=float)
        y_outcome_test = test[outcome_col].to_numpy(dtype=float)

        beta_feature, *_ = np.linalg.lstsq(x_train_conf, y_feature_train, rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train_conf, y_outcome_train, rcond=None)

        resid_feature_train = y_feature_train - x_train_conf @ beta_feature
        resid_outcome_train = y_outcome_train - x_train_conf @ beta_outcome
        resid_feature_test = y_feature_test - x_test_conf @ beta_feature
        resid_outcome_test = y_outcome_test - x_test_conf @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            predicted_resid = float("nan")
            predicted_raw = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            predicted_resid = float(slope * resid_feature_test[0])
            predicted_raw = float((x_test_conf @ beta_outcome)[0] + predicted_resid)

        rows.append(
            {
                id_col: str(test.iloc[0][id_col]),
                outcome_col: float(y_outcome_test[0]),
                "predicted": predicted_raw,
                feature_col: float(y_feature_test[0]),
                "predicted_residual": predicted_resid,
                "actual_residual": float(resid_outcome_test[0]),
            }
        )
        pred_resid.append(predicted_resid)
        act_resid.append(float(resid_outcome_test[0]))

    pred_arr = np.asarray(pred_resid, dtype=float)
    act_arr = np.asarray(act_resid, dtype=float)
    good = np.isfinite(pred_arr) & np.isfinite(act_arr)
    if good.sum() < 3 or np.std(pred_arr[good]) == 0 or np.std(act_arr[good]) == 0:
        corr = float("nan")
    else:
        corr = float(np.corrcoef(pred_arr[good], act_arr[good])[0, 1])

    return corr, pd.DataFrame(rows)


def _evaluate_variation(
    donor_table: pd.DataFrame,
    *,
    variation_name: str,
    n_total: int,
) -> tuple[dict[str, object], pd.DataFrame]:
    feature_col = _feature_column_for(variation_name)
    valid = donor_table.loc[:, ["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUNDS]].dropna().copy()
    analyzable_ids = valid["donor_id"].astype(str).tolist()
    partial = partial_correlation(
        donor_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
    )
    loo_predictive_r, loo_table = _loo_predictions(
        donor_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
        id_col="donor_id",
    )
    influence = leave_one_out_summary(
        donor_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
        id_col="donor_id",
    )

    n_analyzable = int(len(valid))
    partial_r = float(partial.get("partial_r", float("nan")))
    p_value = float(partial.get("p_value", float("nan")))
    selection = _selection_score(partial_r, n_analyzable, n_total)
    gap, penalty, adjusted = _adjusted_metrics(partial_r, loo_predictive_r)

    metrics: dict[str, object] = {
        "variation_name": variation_name,
        "description": VARIATION_BY_NAME[variation_name]["description"],
        "feature_column": feature_col,
        "n_total": int(n_total),
        "n_analyzable": n_analyzable,
        "partial_r": partial_r,
        "p_value": p_value,
        "selection_score": selection,
        "loo_predictive_r": float(loo_predictive_r),
        "is_loo_gap": gap,
        "gap_penalty": penalty,
        "adjusted_score": adjusted,
        "loo_unstable_count": int(influence.get("unstable_count", 0) or 0),
        "loo_max_shift": float(influence.get("max_shift", float("nan"))),
        "donor_ids_used": analyzable_ids,
    }
    return metrics, loo_table


def _sort_value(value: object) -> float:
    try:
        x = float(value)
    except Exception:
        return float("-inf")
    if not math.isfinite(x):
        return float("-inf")
    return x


def _write_report(
    *,
    report_path: Path,
    ranked: list[dict[str, object]],
    best: dict[str, object],
    best_loo_table: pd.DataFrame,
    donor_table: pd.DataFrame,
) -> None:
    winner = str(best["variation_name"])
    loser_rows = [row for row in ranked if row["variation_name"] != winner]
    loser_summary = "; ".join(
        f"{row['variation_name']} selection={float(row['selection_score']):.4f}, partial_r={float(row['partial_r']):.4f}, loo_r={float(row['loo_predictive_r']):.4f}"
        for row in loser_rows
    )
    if not loser_summary:
        loser_summary = "No alternative local variation was tested."

    merged = best_loo_table.merge(
        donor_table.loc[
            :,
            [
                "donor_id",
                "ca1_astrocyte_count",
                "ca1_reactive_astrocyte_count",
                str(best["feature_column"]),
            ],
        ],
        on="donor_id",
        how="left",
    )
    merged["raw_error"] = merged["predicted"] - merged[OUTCOME_COLUMN]
    merged["abs_error"] = merged["raw_error"].abs()
    err_rows = merged.sort_values("abs_error", ascending=False).head(3)

    if len(err_rows) > 0:
        donor_bits = []
        near_threshold = 0
        for row in err_rows.itertuples(index=False):
            astro_n = int(getattr(row, "ca1_astrocyte_count", float("nan")) or 0)
            react_n = int(getattr(row, "ca1_reactive_astrocyte_count", float("nan")) or 0)
            if min(astro_n, react_n) < 2 * MIN_CELLS_PER_GROUP:
                near_threshold += 1
            donor_bits.append(
                f"{row.donor_id} (obs={getattr(row, OUTCOME_COLUMN):.3f}, pred={row.predicted:.3f}, "
                f"feature={getattr(row, str(best['feature_column'])):.3f}, astro={astro_n}, reactive={react_n})"
            )
        error_examples = "; ".join(donor_bits)
        if near_threshold >= 2:
            shared_error = "Two or more of the largest-error donors sit near the minimum-count boundary, so percentile morphology estimates are probably noisier there."
        else:
            shared_error = "The biggest misses are not primarily missing-data cases; they look more like donors where reactive astrocyte size and memory decline are partly decoupled."
    else:
        error_examples = "No analyzable LOO rows were available."
        shared_error = "Error structure could not be assessed."

    if winner == "upper_tail_log_area_shift":
        worked = "The upper-tail version won, which suggests the signal is concentrated in a subset of especially enlarged CA1 reactive astrocytes rather than a uniform shift across all reactive astrocytes."
        failed = "The median-shift version likely diluted that severe-activation subset by averaging over mildly enlarged and more typical reactive astrocytes."
        rationale = "A 75th-percentile shift is biologically coherent for reactive hypertrophy because severe astrocyte activation often appears as a tail of very enlarged cells, not just a small movement of the entire population."
    else:
        worked = "The median-shift version won, which suggests the relevant signal is a broad donor-level enlargement of CA1 reactive astrocytes relative to baseline astrocytes."
        failed = "The upper-tail version appears too tail-sensitive on this cohort size, so a few extreme contours likely added noise instead of sharpening the donor summary."
        rationale = "A median shift is biologically coherent when reactive enlargement is a widespread state change within the CA1 astroglial compartment, making it more stable than a tail-only readout."

    summary = (
        f"Tested the {FAMILY_NAME} family in CA1; {winner} won with selection score "
        f"{float(best['selection_score']):.4f}."
    )

    metrics = f"""## Metrics
Winner `{winner}` (`{best['feature_column']}`):
- IS partial r: {float(best['partial_r']):.4f}
- Selection score: {float(best['selection_score']):.4f}
- LOO predictive r: {float(best['loo_predictive_r']):.4f}
- IS-LOO Gap: {float(best['is_loo_gap']):.4f} (penalty={float(best['gap_penalty']):.4f})
- Adjusted Score: {float(best['adjusted_score']):.4f}
- n_analyzable: {int(best['n_analyzable'])}/{int(best['n_total'])}

Other tested variation ranking:
- {loser_summary}
"""

    findings = f"""## Findings
1. What worked and why: {worked} The donor-level scalar is the within-CA1 log-area shift between Reactive Astrocytes and Astrocytes, so it tracks activation severity rather than raw abundance.
2. What failed and why: {failed}
3. Error pattern: largest held-out errors were {error_examples}. {shared_error}
"""

    rationale_text = f"""## Rationale
{rationale} Restricting the comparison to CA1 sharpens the feature to the same hippocampal field where prior rounds already found pyramidal and reactive-astrocyte signal. Using CA1 Astrocytes as the within-slide reference means the biomarker asks whether reactive astrocytes are unusually enlarged beyond the baseline astrocyte size context of that same donor. That makes it more morphology-specific than abundance or proximity alone. Relative to the current panel, this feature looks biologically adjacent to the accepted CA1 reactive-astrocyte enrichment axis, but it is mechanistically more state-focused, so it still seems plausible as a source of partially nonredundant information.
"""

    interpretation = f"""## Interpretation
Biologically, the signal appears to mean stronger CA1 astrocyte activation marked by reactive-cell hypertrophy. The population is CA1 Reactive Astrocytes referenced against CA1 Astrocytes; the niche is the CA1 field itself; the donor summary is `{winner}` on log contour area; and the simplest observable tissue pattern is a donor having visibly larger reactive astrocyte contours in CA1 than its non-reactive astrocyte baseline.
"""

    next_text = f"""## Next
Stay on the same CA1 reactive-astrocyte morphology axis and test one nearby sweep that changes only the donor aggregation: keep the winning CA1 reactive-vs-astrocyte size contrast, but compare this {winner} statistic with a robust trimmed-upper-tail summary to see whether the winner is capturing a broad state shift or a sparse hypertrophic subpopulation.
"""

    report = "\n".join(
        [
            "## Summary",
            summary,
            "",
            metrics.strip(),
            "",
            findings.strip(),
            "",
            rationale_text.strip(),
            "",
            interpretation.strip(),
            "",
            next_text.strip(),
            "",
        ]
    )
    report_path.write_text(report, encoding="utf-8")


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
    cohort = _cohort_with_covariates(data_root)
    donor_rows = cohort.loc[cohort["donor_id"].astype(str) == str(donor_id)]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    summary = _extract_donor_summary(slide_path)
    variation_name = _canonical_variation(Path(__file__).resolve().parent)
    feature_col = _feature_column_for(variation_name)
    value = summary.get(feature_col, float("nan"))
    if value is None or not np.isfinite(float(value)):
        return None
    return float(value)


def main() -> int:
    base_dir = Path(__file__).resolve().parent
    paths = _paths(base_dir)
    cohort = _cohort_with_covariates("/data").copy()

    donor_summaries: list[dict[str, object]] = []
    for row in cohort.itertuples(index=False):
        slide_path = Path("/data") / str(row.slide_name)
        summary = _extract_donor_summary(slide_path)
        summary["donor_id"] = str(row.donor_id)
        donor_summaries.append(summary)

    summary_df = pd.DataFrame(donor_summaries)
    donor_table = cohort.merge(summary_df, on="donor_id", how="left")

    ranked_rows: list[dict[str, object]] = []
    loo_tables: dict[str, pd.DataFrame] = {}
    for variation in VARIATIONS:
        metrics, loo_table = _evaluate_variation(
            donor_table,
            variation_name=str(variation["name"]),
            n_total=int(len(cohort)),
        )
        ranked_rows.append(metrics)
        loo_tables[str(variation["name"])] = loo_table

    ranked_rows = sorted(
        ranked_rows,
        key=lambda row: (
            _sort_value(row.get("selection_score")),
            _sort_value(abs(float(row.get("partial_r", float("nan"))))),
            _sort_value(row.get("adjusted_score")),
        ),
        reverse=True,
    )
    best = dict(ranked_rows[0])
    best_variation = str(best["variation_name"])
    best_feature_column = str(best["feature_column"])
    best_loo_table = loo_tables[best_variation].copy()

    best_loo_table = best_loo_table.merge(
        donor_table.loc[:, ["donor_id", best_feature_column]],
        on="donor_id",
        how="left",
        suffixes=("", "_full"),
    )
    if f"{best_feature_column}_full" in best_loo_table.columns:
        best_loo_table[best_feature_column] = best_loo_table[f"{best_feature_column}_full"]
        best_loo_table = best_loo_table.drop(columns=[f"{best_feature_column}_full"])

    write_donor_feature_table(
        paths["table"],
        donor_table,
        feature_column=best_feature_column,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUNDS,
        extra_columns=[
            row["feature_column"] for row in VARIATIONS
        ] + ["ca1_astrocyte_count", "ca1_reactive_astrocyte_count", "ca1_total_astroglia_count"],
    )

    state = {
        "best_variation": best_variation,
        "feature_name": FEATURE_NAME,
        "feature_column": best_feature_column,
        "min_cells_per_group": MIN_CELLS_PER_GROUP,
    }
    paths["state"].write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")

    payload = build_results_payload(
        status="ok",
        feature_name=FEATURE_NAME,
        outcome=OUTCOME_COLUMN,
        n_total=int(best["n_total"]),
        n_analyzable=int(best["n_analyzable"]),
        partial_r=float(best["partial_r"]),
        ci_lo=None,
        ci_hi=None,
        p_value=float(best["p_value"]),
        loo_predictive_r=float(best["loo_predictive_r"]),
        loo_unstable_count=int(best["loo_unstable_count"]),
        loo_max_shift=float(best["loo_max_shift"]),
        donor_ids_used=list(best["donor_ids_used"]),
        covariates=CONFOUNDS,
        recomputed_from_raw=True,
        registry_written=False,
        artifacts={
            "donor_feature_table": str(paths["table"]),
            "canonical_state": str(paths["state"]),
            "report": str(paths["report"]),
        },
        best_variation=best_variation,
        feature_column=best_feature_column,
        selection_score=float(best["selection_score"]),
        is_loo_gap=float(best["is_loo_gap"]),
        gap_penalty=float(best["gap_penalty"]),
        adjusted_score=float(best["adjusted_score"]),
        ranked_variations=ranked_rows,
    )
    write_results_payload(paths["results"], payload)

    _write_report(
        report_path=paths["report"],
        ranked=ranked_rows,
        best=best,
        best_loo_table=best_loo_table,
        donor_table=donor_table,
    )

    print(f"HYPOTHESIS FAMILY: {FAMILY_NAME}")
    print(f"BEST VARIATION: {best_variation}")
    print(f"  IS partial r:      {float(best['partial_r']): .4f}")
    print(f"  Selection score:   {float(best['selection_score']): .4f}")
    print(f"  LOO predictive r:  {float(best['loo_predictive_r']): .4f}  (diagnostic)")
    print(f"  IS-LOO Gap:        {float(best['is_loo_gap']): .4f}  (penalty={float(best['gap_penalty']):.4f})")
    print(f"  Adjusted Score:    {float(best['adjusted_score']): .4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for row in ranked_rows:
        print(
            f"  {row['variation_name']}  {float(row['partial_r']): .4f}  "
            f"{float(row['selection_score']): .4f}  {float(row['loo_predictive_r']): .4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    preview = best_loo_table.loc[:, ["donor_id", OUTCOME_COLUMN, "predicted", best_feature_column]].copy()
    for row in preview.itertuples(index=False):
        print(
            f"  {row.donor_id}  {getattr(row, OUTCOME_COLUMN): .4f}  "
            f"{row.predicted: .4f}  {getattr(row, best_feature_column): .4f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
PY
python /scratch/result.py
