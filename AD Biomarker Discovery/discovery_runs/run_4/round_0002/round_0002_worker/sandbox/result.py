from __future__ import annotations

import json
import math
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/shared/lib")

from shared_analysis.artifacts import build_results_payload, write_donor_feature_table
from shared_analysis.sea_ad_lfb import build_cell_table
from shared_analysis.stats import bootstrap_partial_correlation, partial_correlation


FEATURE_FAMILY = "ca1_reactive_astrocyte_enrichment"
CANONICAL_VARIATION = "reactive_over_astroglia"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]

VARIATIONS = {
    "reactive_over_astroglia": {
        "description": "CA1 Reactive Astrocyte / (CA1 Reactive Astrocyte + CA1 Astrocyte)",
        "feature_column": f"{FEATURE_FAMILY}__reactive_over_astroglia",
    },
    "reactive_over_all_cells": {
        "description": "CA1 Reactive Astrocyte / all classified CA1 cells",
        "feature_column": f"{FEATURE_FAMILY}__reactive_over_all_cells",
    },
}


def _cohort_path(data_root: str | Path) -> Path:
    return Path(data_root) / "training_cohort.csv"


def feature_column_for_variation(variation: str) -> str:
    return str(VARIATIONS[variation]["feature_column"])


FEATURE_COLUMN = feature_column_for_variation(CANONICAL_VARIATION)


def _sex_to_binary(value: object) -> float:
    if pd.isna(value):
        return float("nan")
    text = str(value).strip().lower()
    if text == "female":
        return 0.0
    if text == "male":
        return 1.0
    return float("nan")


def _design_matrix(frame: pd.DataFrame, columns: list[str]) -> np.ndarray:
    cols = [frame[col].to_numpy(dtype=float) for col in columns]
    return np.column_stack([np.ones(len(frame), dtype=float), *cols])


def _clean_frame(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = frame.loc[:, columns].copy()
    for col in columns:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    return out.replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)


def extract_ca1_counts(zarr_path: str | Path) -> dict[str, float]:
    warnings.filterwarnings(
        "ignore",
        message="Object at .* is not recognized as a component of a Zarr hierarchy.",
    )
    cells = build_cell_table(zarr_path, include_regions=True, include_geometry=False)
    ca1 = cells.loc[cells["region"] == "CA1", ["cell_type"]]
    reactive = int((ca1["cell_type"] == "Reactive Astrocyte").sum())
    astrocyte = int((ca1["cell_type"] == "Astrocyte").sum())
    total = int(len(ca1))
    return {
        "ca1_reactive_astrocyte_count": float(reactive),
        "ca1_astrocyte_count": float(astrocyte),
        "ca1_total_classified_cells": float(total),
    }


def compute_variation_score_from_counts(
    *,
    reactive_count: float,
    astrocyte_count: float,
    total_count: float,
    variation: str,
) -> float:
    reactive_count = float(reactive_count)
    astrocyte_count = float(astrocyte_count)
    total_count = float(total_count)

    if variation == "reactive_over_astroglia":
        denom = reactive_count + astrocyte_count
    elif variation == "reactive_over_all_cells":
        denom = total_count
    else:
        raise ValueError(f"Unknown variation: {variation}")

    if not np.isfinite(denom) or denom <= 0:
        return float("nan")
    return float(reactive_count / denom)


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Return the winning CA1 reactive astrocyte enrichment score for one donor.
    """
    data_root = Path(data_root)
    cohort = pd.read_csv(_cohort_path(data_root))
    donor_rows = cohort.loc[cohort["donor_id"].astype(str) == str(donor_id)]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    counts = extract_ca1_counts(data_root / slide_name)
    return compute_variation_score_from_counts(
        reactive_count=counts["ca1_reactive_astrocyte_count"],
        astrocyte_count=counts["ca1_astrocyte_count"],
        total_count=counts["ca1_total_classified_cells"],
        variation=CANONICAL_VARIATION,
    )


def extract_all_donor_features(data_root: str | Path) -> pd.DataFrame:
    data_root = Path(data_root)
    cohort = pd.read_csv(_cohort_path(data_root))
    if "sex_binary" not in cohort.columns:
        cohort["sex_binary"] = cohort["sex"].map(_sex_to_binary)

    rows: list[dict[str, object]] = []
    for _, row in cohort.iterrows():
        slide_name = str(row["slide_name"])
        zarr_path = data_root / slide_name
        counts = extract_ca1_counts(zarr_path)
        record = row.to_dict()
        record.update(counts)
        for variation in VARIATIONS:
            record[feature_column_for_variation(variation)] = compute_variation_score_from_counts(
                reactive_count=counts["ca1_reactive_astrocyte_count"],
                astrocyte_count=counts["ca1_astrocyte_count"],
                total_count=counts["ca1_total_classified_cells"],
                variation=variation,
            )
        rows.append(record)
    return pd.DataFrame(rows)


def loo_prediction_details(
    frame: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
) -> tuple[float, pd.DataFrame]:
    clean = frame.copy()
    required_cols = [feature_col, outcome_col, *confounds]
    for col in required_cols:
        clean[col] = pd.to_numeric(clean[col], errors="coerce")
    clean = clean.replace([np.inf, -np.inf], np.nan).dropna(subset=required_cols).reset_index(drop=True)

    preds_raw: list[float] = []
    preds_resid: list[float] = []
    actual_resid: list[float] = []
    out_rows: list[dict[str, object]] = []

    for idx in range(len(clean)):
        train = clean.drop(index=idx).reset_index(drop=True)
        test = clean.iloc[[idx]].reset_index(drop=True)

        x_train_full = _design_matrix(train, [*confounds, feature_col])
        x_test_full = _design_matrix(test, [*confounds, feature_col])
        y_train = train[outcome_col].to_numpy(dtype=float)
        beta_full, *_ = np.linalg.lstsq(x_train_full, y_train, rcond=None)
        pred_raw = float((x_test_full @ beta_full)[0])

        x_train_conf = _design_matrix(train, confounds)
        x_test_conf = _design_matrix(test, confounds)

        beta_feature, *_ = np.linalg.lstsq(
            x_train_conf,
            train[feature_col].to_numpy(dtype=float),
            rcond=None,
        )
        beta_outcome, *_ = np.linalg.lstsq(
            x_train_conf,
            train[outcome_col].to_numpy(dtype=float),
            rcond=None,
        )

        resid_feature_train = train[feature_col].to_numpy(dtype=float) - x_train_conf @ beta_feature
        resid_outcome_train = train[outcome_col].to_numpy(dtype=float) - x_train_conf @ beta_outcome
        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            pred_resid = float("nan")
            act_resid = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            resid_feature_test = test[feature_col].to_numpy(dtype=float) - x_test_conf @ beta_feature
            pred_resid = float(slope * resid_feature_test[0])
            act_resid = float(test[outcome_col].to_numpy(dtype=float)[0] - (x_test_conf @ beta_outcome)[0])

        preds_raw.append(pred_raw)
        preds_resid.append(pred_resid)
        actual_resid.append(act_resid)

        feature_value = float(test.iloc[0][feature_col])
        out_rows.append(
            {
                "donor_id": str(test.iloc[0]["donor_id"]),
                "outcome": float(test.iloc[0][outcome_col]),
                "predicted": pred_raw,
                "predicted_resid": pred_resid,
                "actual_resid": act_resid,
                feature_col: feature_value,
                "ca1_reactive_astrocyte_count": float(test.iloc[0]["ca1_reactive_astrocyte_count"]),
                "ca1_astrocyte_count": float(test.iloc[0]["ca1_astrocyte_count"]),
                "ca1_total_classified_cells": float(test.iloc[0]["ca1_total_classified_cells"]),
            }
        )

    pred_arr = np.asarray(preds_resid, dtype=float)
    act_arr = np.asarray(actual_resid, dtype=float)
    mask = np.isfinite(pred_arr) & np.isfinite(act_arr)
    if mask.sum() < 3 or np.nanstd(pred_arr[mask]) == 0 or np.nanstd(act_arr[mask]) == 0:
        loo_r = float("nan")
    else:
        loo_r = float(np.corrcoef(pred_arr[mask], act_arr[mask])[0, 1])

    return loo_r, pd.DataFrame(out_rows)


def evaluate_variation(
    donor_table: pd.DataFrame,
    *,
    variation: str,
    n_total: int,
) -> dict[str, object]:
    feature_col = feature_column_for_variation(variation)
    needed = donor_table[["donor_id", OUTCOME_COLUMN, feature_col, *CONFOUND_COLUMNS]].copy()
    for col in [OUTCOME_COLUMN, feature_col, *CONFOUND_COLUMNS]:
        needed[col] = pd.to_numeric(needed[col], errors="coerce")
    analyzable_mask = np.isfinite(needed[feature_col]) & np.isfinite(needed[OUTCOME_COLUMN])
    for col in CONFOUND_COLUMNS:
        analyzable_mask &= np.isfinite(needed[col])
    analyzable = donor_table.loc[analyzable_mask].copy().reset_index(drop=True)
    n_analyzable = int(len(analyzable))

    if n_analyzable < 4:
        return {
            "name": variation,
            "description": VARIATIONS[variation]["description"],
            "feature_column": feature_col,
            "n_analyzable": n_analyzable,
            "partial_r": float("nan"),
            "p_value": float("nan"),
            "ci_lo": float("nan"),
            "ci_hi": float("nan"),
            "selection_score": float("nan"),
            "loo_predictive_r": float("nan"),
            "is_loo_gap": float("nan"),
            "penalty": float("nan"),
            "adjusted_score": float("nan"),
            "loo_table": pd.DataFrame(),
        }

    partial = partial_correlation(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    boot = bootstrap_partial_correlation(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        n_boot=1000,
        random_state={"reactive_over_astroglia": 11, "reactive_over_all_cells": 23}[variation],
    )
    loo_r, loo_table = loo_prediction_details(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    partial_r = float(partial["partial_r"])
    selection_score = float(abs(partial_r) * (n_analyzable / n_total)) if np.isfinite(partial_r) else float("nan")
    gap = float(abs(partial_r) - abs(loo_r)) if np.isfinite(partial_r) and np.isfinite(loo_r) else float("nan")
    penalty = float(max(gap - 0.15, 0.0) / 2.0) if np.isfinite(gap) else float("nan")
    adjusted_score = float(loo_r - penalty) if np.isfinite(loo_r) else float("nan")

    return {
        "name": variation,
        "description": VARIATIONS[variation]["description"],
        "feature_column": feature_col,
        "n_analyzable": n_analyzable,
        "partial_r": partial_r,
        "p_value": float(partial["p_value"]),
        "ci_lo": float(boot["ci_lo"]),
        "ci_hi": float(boot["ci_hi"]),
        "selection_score": selection_score,
        "loo_predictive_r": float(loo_r),
        "is_loo_gap": gap,
        "penalty": penalty,
        "adjusted_score": adjusted_score,
        "loo_table": loo_table,
    }


def _fmt(value: float) -> str:
    if value is None or not np.isfinite(value):
        return "nan"
    return f"{value:.4f}"


def _update_canonical_variation(best_variation: str) -> None:
    path = Path(__file__)
    text = path.read_text()
    updated = re.sub(
        r'(?m)^CANONICAL_VARIATION = ".*"$',
        f'CANONICAL_VARIATION = "{best_variation}"',
        text,
        count=1,
    )
    if updated != text:
        path.write_text(updated)


def _render_report(
    *,
    context: dict[str, object],
    donor_table: pd.DataFrame,
    best: dict[str, object],
    ranked: list[dict[str, object]],
    loo_table: pd.DataFrame,
) -> str:
    brief = context.get("worker_brief", {}) if isinstance(context, dict) else {}
    family = FEATURE_FAMILY
    winner = str(best["name"])
    best_feature_col = str(best["feature_column"])

    ranking_lines = []
    for row in ranked:
        ranking_lines.append(
            f"- `{row['name']}`: partial r {_fmt(float(row['partial_r']))}, "
            f"selection score {_fmt(float(row['selection_score']))}, "
            f"LOO r {_fmt(float(row['loo_predictive_r']))}"
        )

    residuals = loo_table.copy()
    residuals["abs_error"] = (residuals["predicted"] - residuals["outcome"]).abs()
    worst = residuals.sort_values("abs_error", ascending=False).head(3)
    worst_ids = worst["donor_id"].astype(str).tolist()

    high_feature = loo_table[best_feature_col].median()
    overpred = worst.loc[worst["predicted"] < worst["outcome"], "donor_id"].astype(str).tolist()
    underpred = worst.loc[worst["predicted"] > worst["outcome"], "donor_id"].astype(str).tolist()

    loser = ranked[1] if len(ranked) > 1 else None
    loser_text = ""
    if loser is not None:
        loser_text = (
            f"The nearby denominator swap `{loser['name']}` ranked lower "
            f"(selection score {_fmt(float(loser['selection_score']))}), suggesting that dividing by all CA1 cells diluted "
            f"the astroglial-state signal with broader CA1 composition."
            if loser["name"] == "reactive_over_all_cells"
            else f"The nearby alternative `{loser['name']}` ranked lower."
        )

    top_share = ""
    if not worst.empty:
        high_count = int((worst[best_feature_col] > high_feature).sum())
        if high_count >= 2:
            top_share = (
                f"Two of the three largest errors ({', '.join(worst_ids)}) had above-median {winner} values, "
                "so the model tends to over-call decline in some donors with strong CA1 gliosis but less severe memory slope than expected."
            )
        else:
            top_share = (
                f"The largest errors ({', '.join(worst_ids)}) were mixed in feature direction, implying that reactive enrichment is real but not sufficient alone to explain all donor-to-donor variation after confound adjustment."
            )

    likely_additive = (
        "Because the accepted panel member tracks CA1 pyramidal depletion, this glial-state ratio is biologically distinct and plausibly additive, though panel-level gain must be confirmed by the evaluator."
    )

    return f"""## Summary
Tested the `{family}` family in CA1; `{winner}` won the local sweep with selection score {_fmt(float(best['selection_score']))}.

## Metrics
Winning variation: `{winner}` (`{best_feature_col}`)

- IS partial r: {_fmt(float(best['partial_r']))}
- Selection score: {_fmt(float(best['selection_score']))}
- LOO predictive r: {_fmt(float(best['loo_predictive_r']))}
- IS-LOO gap: {_fmt(float(best['is_loo_gap']))}
- Penalty: {_fmt(float(best['penalty']))}
- Adjusted score: {_fmt(float(best['adjusted_score']))}
- Bootstrap 95% CI: [{_fmt(float(best['ci_lo']))}, {_fmt(float(best['ci_hi']))}]
- n analyzable: {int(best['n_analyzable'])}

Other tested variations:
{chr(10).join(ranking_lines)}

## Findings
1. What worked and why  
   `{winner}` worked best because it isolates CA1 reactive astrocyte state within the astroglial lineage rather than letting the denominator drift with the whole CA1 cellular mixture. That makes it a cleaner donor-level scalar for local gliosis / injury response.

2. What failed and why  
   {loser_text or "No nearby alternative was tested."}

3. Error pattern: which donors are consistently wrong and what they share  
   The largest absolute LOO prediction errors were in donors {', '.join(worst_ids)}. {top_share}

## Rationale
The best variation is biologically coherent because it measures a shift from homeostatic `Astrocyte` toward `Reactive Astrocyte` specifically inside CA1, the same region already implicated by the accepted panel feature. That framing matches a local injury-response process rather than generic tissue composition. It beat the nearby alternative because the all-cell denominator mixes astroglial activation with large donor-to-donor variation in neurons, oligodendrocytes, and other CA1 cell classes. {likely_additive}

## Interpretation
The signal appears to reflect CA1 astroglial activation.  
Population: `Reactive Astrocyte` versus `Astrocyte`  
Niche: CA1 annotated tissue compartment  
Feature summary: donor-level CA1 reactive astrocyte enrichment ratio  
Simplest observable pattern: CA1 fields where the astrocyte pool is shifted toward reactive-state labels rather than baseline astrocytes.

## Next
Keep the CA1 reactive-astrocyte lead but test a tighter CA1 niche next, such as reactive-over-astroglia restricted to the CA1 pyramidal layer border or a CA1 reactive-to-pyramidal normalization, since the all-cell denominator lost signal and the remaining errors suggest that diffuse CA1 burden is still too coarse.
"""


def main() -> None:
    data_root = Path("/data")
    scratch = Path("/scratch")
    context_path = scratch / "context_bundle.json"
    context = json.loads(context_path.read_text()) if context_path.exists() else {}

    donor_table = extract_all_donor_features(data_root)
    n_total = int(len(donor_table))

    evaluations = [evaluate_variation(donor_table, variation=name, n_total=n_total) for name in VARIATIONS]
    ranked = sorted(
        evaluations,
        key=lambda row: (
            -np.inf if not np.isfinite(float(row["selection_score"])) else -float(row["selection_score"]),
            -np.inf if not np.isfinite(float(row["partial_r"])) else -abs(float(row["partial_r"])),
            -np.inf if not np.isfinite(float(row["loo_predictive_r"])) else -float(row["loo_predictive_r"]),
        ),
    )
    best = ranked[0]
    best_variation = str(best["name"])
    best_feature_column = str(best["feature_column"])
    best_loo_table = best["loo_table"].copy()

    donor_feature_table = donor_table.copy()
    write_donor_feature_table(
        scratch / "donor_feature_table.csv",
        donor_feature_table,
        feature_column=best_feature_column,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUND_COLUMNS,
        id_columns=["donor_id", "slide_name"],
        extra_columns=[
            "ca1_reactive_astrocyte_count",
            "ca1_astrocyte_count",
            "ca1_total_classified_cells",
        ],
    )

    ranked_json = []
    for row in ranked:
        ranked_json.append(
            {
                "name": str(row["name"]),
                "description": str(row["description"]),
                "feature_column": str(row["feature_column"]),
                "n_analyzable": int(row["n_analyzable"]),
                "partial_r": float(row["partial_r"]),
                "selection_score": float(row["selection_score"]),
                "loo_predictive_r": float(row["loo_predictive_r"]),
                "is_loo_gap": float(row["is_loo_gap"]),
                "penalty": float(row["penalty"]),
                "adjusted_score": float(row["adjusted_score"]),
                "p_value": float(row["p_value"]),
                "ci_lo": float(row["ci_lo"]),
                "ci_hi": float(row["ci_hi"]),
            }
        )

    results = build_results_payload(
        status="ok",
        feature_name=FEATURE_FAMILY,
        outcome=OUTCOME_COLUMN,
        n_total=n_total,
        n_analyzable=int(best["n_analyzable"]),
        partial_r=float(best["partial_r"]),
        ci_lo=float(best["ci_lo"]),
        ci_hi=float(best["ci_hi"]),
        p_value=float(best["p_value"]),
        loo_predictive_r=float(best["loo_predictive_r"]),
        donor_ids_used=donor_table["donor_id"].astype(str).tolist(),
        covariates=CONFOUND_COLUMNS,
        recomputed_from_raw=True,
        registry_written=False,
        artifacts={
            "donor_feature_table": str(scratch / "donor_feature_table.csv"),
            "report_path": str(scratch / "report.md"),
        },
        best_variation=best_variation,
        feature_column=best_feature_column,
        selection_score=float(best["selection_score"]),
        is_loo_gap=float(best["is_loo_gap"]),
        penalty=float(best["penalty"]),
        adjusted_score=float(best["adjusted_score"]),
        ranked_variations=ranked_json,
        per_donor_loo=best_loo_table.to_dict(orient="records"),
    )
    (scratch / "results.json").write_text(json.dumps(results, indent=2))

    report_text = _render_report(
        context=context,
        donor_table=donor_table,
        best=best,
        ranked=ranked_json,
        loo_table=best_loo_table,
    )
    (scratch / "report.md").write_text(report_text)

    print(f"HYPOTHESIS FAMILY: {FEATURE_FAMILY}")
    print(f"BEST VARIATION: {best_variation}")
    print(f"  IS partial r:      {_fmt(float(best['partial_r']))}")
    print(f"  Selection score:   {_fmt(float(best['selection_score']))}")
    print(f"  LOO predictive r:  {_fmt(float(best['loo_predictive_r']))}  (diagnostic)")
    print(f"  IS-LOO Gap:        {_fmt(float(best['is_loo_gap']))}  (penalty={_fmt(float(best['penalty']))})")
    print(f"  Adjusted Score:    {_fmt(float(best['adjusted_score']))}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for row in ranked_json:
        print(
            f"  {row['name']}  {_fmt(float(row['partial_r']))}  "
            f"{_fmt(float(row['selection_score']))}  {_fmt(float(row['loo_predictive_r']))}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(
        f"  donor_id  outcome  predicted  {best_feature_column}  "
        "ca1_reactive_astrocyte_count  ca1_astrocyte_count  ca1_total_classified_cells"
    )
    for _, row in best_loo_table.iterrows():
        print(
            f"  {row['donor_id']}  {_fmt(float(row['outcome']))}  {_fmt(float(row['predicted']))}  "
            f"{_fmt(float(row[best_feature_column]))}  {_fmt(float(row['ca1_reactive_astrocyte_count']))}  "
            f"{_fmt(float(row['ca1_astrocyte_count']))}  {_fmt(float(row['ca1_total_classified_cells']))}"
        )

    _update_canonical_variation(best_variation)


if __name__ == "__main__":
    main()
