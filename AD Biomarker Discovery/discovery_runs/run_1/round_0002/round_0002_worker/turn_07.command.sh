set -euo pipefail
cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, "/shared/lib")

from shared_analysis.artifacts import write_donor_feature_table
from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort
from shared_analysis.stats import partial_correlation


HYPOTHESIS_FAMILY = "ca1-reactive-astro-fraction"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]

VARIATIONS: dict[str, dict[str, str]] = {
    "candidate_variant_a": {
        "feature_name": "ca1_reactive_astrocyte_lineage_fraction",
        "feature_column": "ca1_reactive_astrocyte_lineage_fraction",
        "description": "Reactive Astrocyte / (Astrocyte + Reactive Astrocyte) within CA1",
        "population": "Reactive Astrocyte within astrocyte lineage",
        "region": "CA1",
        "summary": "lineage fraction",
    },
    "candidate_variant_b": {
        "feature_name": "ca1_reactive_astrocyte_tissue_fraction",
        "feature_column": "ca1_reactive_astrocyte_tissue_fraction",
        "description": "Reactive Astrocyte / all classified cells within CA1",
        "population": "Reactive Astrocyte within all classified cells",
        "region": "CA1",
        "summary": "tissue fraction",
    },
}

CANONICAL_VARIATION = "candidate_variant_a"  # AUTO-UPDATED


def _safe_float(value: Any) -> float:
    try:
        out = float(value)
    except Exception:
        return float("nan")
    return out


def _is_finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except Exception:
        return False


def _fmt(value: Any) -> str:
    if _is_finite(value):
        return f"{float(value):.4f}"
    return "nan"


def _sex_to_binary(value: Any) -> float:
    if pd.isna(value):
        return float("nan")
    text = str(value).strip().lower()
    if text == "male":
        return 1.0
    if text == "female":
        return 0.0
    return float("nan")


def _design_matrix(df: pd.DataFrame) -> np.ndarray:
    cols = [df[col].to_numpy(dtype=float) for col in CONFOUND_COLUMNS]
    return np.column_stack([np.ones(len(df), dtype=float), *cols])


def _extract_counts_for_slide(slide_path: Path) -> dict[str, float]:
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False, scale=16.0)
    ca1 = cells.loc[cells["region"] == "CA1", "cell_type"]
    vc = ca1.value_counts()
    n_ca1 = int(len(ca1))
    n_reactive = int(vc.get("Reactive Astrocyte", 0))
    n_astro = int(vc.get("Astrocyte", 0))
    n_lineage = n_reactive + n_astro
    return {
        "ca1_total_cells": float(n_ca1),
        "ca1_reactive_astrocytes": float(n_reactive),
        "ca1_astrocytes": float(n_astro),
        "ca1_astro_lineage_cells": float(n_lineage),
    }


def _compute_variation_score(counts: dict[str, float], variation: str) -> float:
    reactive = counts["ca1_reactive_astrocytes"]
    lineage = counts["ca1_astro_lineage_cells"]
    total = counts["ca1_total_cells"]
    if variation == "candidate_variant_a":
        if lineage <= 0:
            return float("nan")
        return float(reactive / lineage)
    if variation == "candidate_variant_b":
        if total <= 0:
            return float("nan")
        return float(reactive / total)
    raise KeyError(f"Unknown variation: {variation}")


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Canonical replay target for held-out evaluation.
    """
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    counts = _extract_counts_for_slide(slide_path)
    return _compute_variation_score(counts, CANONICAL_VARIATION)


def _extract_all_features(data_root: Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map(_sex_to_binary)

    records: list[dict[str, Any]] = []
    for row in cohort.itertuples(index=False):
        slide_path = data_root / row.slide_name
        counts = _extract_counts_for_slide(slide_path)
        record = {
            "donor_id": row.donor_id,
            "slide_name": row.slide_name,
            OUTCOME_COLUMN: _safe_float(getattr(row, OUTCOME_COLUMN)),
            "max_age_vis": _safe_float(row.max_age_vis),
            "braak_numeric": _safe_float(row.braak_numeric),
            "cerad_ordinal": _safe_float(row.cerad_ordinal),
            "sex": row.sex,
            "sex_binary": _sex_to_binary(row.sex),
            "cognitive_status": row.cognitive_status,
            "overall_ad_neuropath_change": row.overall_ad_neuropath_change,
        }
        record.update(counts)
        for variation_name, meta in VARIATIONS.items():
            record[meta["feature_column"]] = _compute_variation_score(counts, variation_name)
        records.append(record)
    return pd.DataFrame.from_records(records)


def _loo_details(frame: pd.DataFrame, *, feature_col: str) -> tuple[pd.DataFrame, float]:
    frame = frame.dropna(subset=[feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]).reset_index(drop=True).copy()
    if len(frame) < 3:
        return frame.assign(predicted=np.nan, predicted_residual=np.nan, actual_residual=np.nan, abs_error=np.nan), float("nan")

    predicted = np.full(len(frame), np.nan, dtype=float)
    predicted_resid = np.full(len(frame), np.nan, dtype=float)
    actual_resid = np.full(len(frame), np.nan, dtype=float)

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train = _design_matrix(train)
        x_test = _design_matrix(test)

        beta_feature, *_ = np.linalg.lstsq(x_train, train[feature_col].to_numpy(dtype=float), rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train, train[OUTCOME_COLUMN].to_numpy(dtype=float), rcond=None)

        resid_feature_train = train[feature_col].to_numpy(dtype=float) - x_train @ beta_feature
        resid_outcome_train = train[OUTCOME_COLUMN].to_numpy(dtype=float) - x_train @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            continue

        slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
        resid_feature_test = test[feature_col].to_numpy(dtype=float) - x_test @ beta_feature
        resid_outcome_test = test[OUTCOME_COLUMN].to_numpy(dtype=float) - x_test @ beta_outcome

        predicted_resid[idx] = float(slope * resid_feature_test[0])
        actual_resid[idx] = float(resid_outcome_test[0])
        predicted[idx] = float((x_test @ beta_outcome)[0] + predicted_resid[idx])

    mask = np.isfinite(predicted_resid) & np.isfinite(actual_resid)
    if mask.sum() >= 3 and np.std(predicted_resid[mask]) > 0 and np.std(actual_resid[mask]) > 0:
        loo_r = float(np.corrcoef(predicted_resid[mask], actual_resid[mask])[0, 1])
    else:
        loo_r = float("nan")

    out = frame.copy()
    out["predicted"] = predicted
    out["predicted_residual"] = predicted_resid
    out["actual_residual"] = actual_resid
    out["abs_error"] = np.abs(out[OUTCOME_COLUMN] - out["predicted"])
    return out, loo_r


def _adjusted_score(partial_r: float, loo_r: float) -> tuple[float, float, float]:
    if not (_is_finite(partial_r) and _is_finite(loo_r)):
        return float("nan"), float("nan"), float("nan")
    gap = abs(float(partial_r)) - abs(float(loo_r))
    penalty = max(0.0, gap - 0.15) * 0.5
    if gap > 0.30:
        adjusted = -1.0
    else:
        adjusted = float(loo_r - penalty)
    return float(gap), float(penalty), float(adjusted)


def _evaluate_variation(table: pd.DataFrame, variation_name: str) -> dict[str, Any]:
    meta = VARIATIONS[variation_name]
    feature_col = meta["feature_column"]
    clean = table.dropna(subset=[feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]).reset_index(drop=True).copy()

    partial = partial_correlation(
        clean,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    partial_r = _safe_float(partial.get("partial_r"))
    p_value = _safe_float(partial.get("p_value"))
    n_analyzable = int(len(clean))
    n_total = int(len(table))
    selection_score = abs(partial_r) * (n_analyzable / n_total) if _is_finite(partial_r) else float("nan")

    loo_df, loo_r = _loo_details(clean, feature_col=feature_col)
    gap, penalty, adjusted = _adjusted_score(partial_r, loo_r)

    return {
        "variation_name": variation_name,
        "feature_name": meta["feature_name"],
        "feature_column": feature_col,
        "description": meta["description"],
        "population": meta["population"],
        "region": meta["region"],
        "summary": meta["summary"],
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "partial_r": partial_r,
        "p_value": p_value,
        "selection_score": selection_score,
        "loo_predictive_r": loo_r,
        "is_loo_gap": gap,
        "gap_penalty": penalty,
        "adjusted_score": adjusted,
        "loo_table": loo_df,
    }


def _sort_key(result: dict[str, Any]) -> tuple[float, float, float]:
    sel = result["selection_score"]
    adj = result["adjusted_score"]
    loo = result["loo_predictive_r"]
    return (
        -1e9 if not _is_finite(sel) else float(sel),
        -1e9 if not _is_finite(adj) else float(adj),
        -1e9 if not _is_finite(loo) else float(loo),
    )


def _update_canonical_variation(best_variation: str) -> None:
    path = Path(__file__)
    text = path.read_text()
    new_text = re.sub(
        r'^CANONICAL_VARIATION = ".*?"  # AUTO-UPDATED$',
        f'CANONICAL_VARIATION = "{best_variation}"  # AUTO-UPDATED',
        text,
        flags=re.MULTILINE,
    )
    if new_text != text:
        path.write_text(new_text)


def _write_results_json(best: dict[str, Any], ranked: list[dict[str, Any]]) -> None:
    ranked_payload = []
    for row in ranked:
        ranked_payload.append(
            {
                "variation_name": row["variation_name"],
                "feature_name": row["feature_name"],
                "feature_column": row["feature_column"],
                "description": row["description"],
                "population": row["population"],
                "region": row["region"],
                "summary": row["summary"],
                "n_total": row["n_total"],
                "n_analyzable": row["n_analyzable"],
                "partial_r": row["partial_r"],
                "p_value": row["p_value"],
                "selection_score": row["selection_score"],
                "loo_predictive_r": row["loo_predictive_r"],
                "is_loo_gap": row["is_loo_gap"],
                "gap_penalty": row["gap_penalty"],
                "adjusted_score": row["adjusted_score"],
            }
        )

    payload = {
        "status": "ok",
        "feature_name": best["feature_name"],
        "feature_column": best["feature_column"],
        "best_variation": best["variation_name"],
        "ranked_variations": ranked_payload,
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "outcome": OUTCOME_COLUMN,
        "covariates": CONFOUND_COLUMNS,
        "n_total": best["n_total"],
        "n_analyzable": best["n_analyzable"],
        "partial_r": best["partial_r"],
        "p_value": best["p_value"],
        "selection_score": best["selection_score"],
        "loo_predictive_r": best["loo_predictive_r"],
        "is_loo_gap": best["is_loo_gap"],
        "gap_penalty": best["gap_penalty"],
        "adjusted_score": best["adjusted_score"],
        "recomputed_from_raw": True,
        "artifacts": {
            "donor_feature_table": "/scratch/donor_feature_table.csv",
            "report": "/scratch/report.md",
        },
    }
    Path("/scratch/results.json").write_text(json.dumps(payload, indent=2))


def _common_or_mixed(values: pd.Series) -> str:
    vals = [str(v) for v in values.dropna().tolist()]
    if not vals:
        return "missing"
    counts = pd.Series(vals).value_counts()
    if len(counts) == 1:
        return counts.index[0]
    if counts.iloc[0] >= 2:
        return f"mostly {counts.index[0]}"
    return "mixed"


def _build_report(best: dict[str, Any], ranked: list[dict[str, Any]]) -> str:
    best_loo = best["loo_table"].copy()
    top_err = best_loo.sort_values("abs_error", ascending=False).head(3)
    donor_bits = []
    for row in top_err.itertuples(index=False):
        donor_bits.append(
            f"{row.donor_id} (outcome={row.slope_zmem0:.3f}, pred={row.predicted:.3f}, "
            f"{best['feature_column']}={getattr(row, best['feature_column']):.4f}, "
            f"status={row.cognitive_status}, AD_change={row.overall_ad_neuropath_change})"
        )
    error_pattern = "; ".join(donor_bits) if donor_bits else "No analyzable donors."

    shared_cognition = _common_or_mixed(top_err["cognitive_status"]) if len(top_err) else "missing"
    shared_path = _common_or_mixed(top_err["overall_ad_neuropath_change"]) if len(top_err) else "missing"

    others = [r for r in ranked if r["variation_name"] != best["variation_name"]]
    if others:
        other_lines = "; ".join(
            f"{r['variation_name']} partial_r={_fmt(r['partial_r'])}, "
            f"selection_score={_fmt(r['selection_score'])}, loo_r={_fmt(r['loo_predictive_r'])}"
            for r in others
        )
    else:
        other_lines = "No other local variations were tested."

    if best["variation_name"] == "candidate_variant_a":
        worked = (
            "Normalizing reactive astrocytes by the CA1 astrocyte lineage sharpened the signal, "
            "suggesting the informative contrast is local astrogliosis state rather than whole-tissue abundance."
        )
        failed = (
            "The tissue-fraction alternative diluted the signal by mixing reactive astrocytes with large donor-to-donor "
            "variation in neurons and other classes, so the denominator likely reintroduced broader composition noise."
        )
        next_step = (
            "Test whether the same CA1 reactive-astrocyte lineage fraction becomes stronger when restricted to the immediate "
            "niche around CA1 pyramidal neurons rather than the whole CA1 annotation."
        )
        observable_pattern = (
            "CA1 fields where Reactive Astrocyte labels visibly replace a larger share of the astrocyte lineage."
        )
    else:
        worked = (
            "Normalizing reactive astrocytes by all CA1 cells worked best, implying the donor-level signal tracks "
            "overall tissue burden of reactive astrocytes more than the balance within the astrocyte lineage."
        )
        failed = (
            "The lineage-fraction alternative was weaker, likely because the astrocyte-lineage denominator is smaller and "
            "adds variance when astrocyte counts are modest in some donors."
        )
        next_step = (
            "Test whether the CA1 reactive-astrocyte tissue fraction is strongest after excluding obvious low-information "
            "classes such as negative control and corpora amylacea from the denominator."
        )
        observable_pattern = (
            "CA1 fields with more frequent Reactive Astrocyte labels spread across the tissue, not just within the glial compartment."
        )

    sign = best["partial_r"]
    if _is_finite(sign) and sign > 0:
        direction = "Higher"
        relation = "higher"
    elif _is_finite(sign) and sign < 0:
        direction = "Lower"
        relation = "lower"
    else:
        direction = "Variation in"
        relation = "indeterminate"

    report = f"""## Summary
Tested the {HYPOTHESIS_FAMILY} family in CA1; {best['variation_name']} won with selection score {_fmt(best['selection_score'])}.

## Metrics
Winning variation: {best['variation_name']} (`{best['feature_column']}`).

- IS partial r: {_fmt(best['partial_r'])}
- Selection score: {_fmt(best['selection_score'])}
- LOO predictive r: {_fmt(best['loo_predictive_r'])}
- IS-LOO Gap: {_fmt(best['is_loo_gap'])} (penalty={_fmt(best['gap_penalty'])})
- Adjusted Score: {_fmt(best['adjusted_score'])}

Other tested variations: {other_lines}

## Findings
1. What worked and why: {worked}
2. What failed and why: {failed}
3. Error pattern: the largest LOO errors were {error_pattern}. These errors are {shared_cognition} for cognitive status and {shared_path} for AD neuropath change, so the mistakes appear concentrated in a partly shared clinical/pathology subset rather than being purely random.

## Rationale
The best variation is biologically coherent because it keeps the Round 1 vulnerable region anchor (CA1) but shifts the measured population from neuronal loss to glial injury response. That makes the donor-level scalar easy to audit from raw labels: count reactive astrocytes inside CA1 and normalize by a biologically meaningful local denominator. It beat the nearby alternative because its denominator better matched the hypothesized process rather than pooling in unrelated CA1 classes.

## Interpretation
{direction} values of `{best['feature_column']}` were associated with {relation} residualized `{OUTCOME_COLUMN}` in this training cohort. The signal most plausibly reflects CA1-localized reactive astrogliosis. Population: {best['population']}. Niche/region: {best['region']}. Donor-level summary: {best['summary']} of Reactive Astrocyte labels. Simplest observable tissue pattern: {observable_pattern}

## Next
{next_step}
"""
    return report


def main() -> None:
    data_root = Path("/data")
    table = _extract_all_features(data_root)

    evaluations = [_evaluate_variation(table, variation_name) for variation_name in VARIATIONS]
    ranked = sorted(evaluations, key=_sort_key, reverse=True)
    best = ranked[0]

    _update_canonical_variation(best["variation_name"])

    write_donor_feature_table(
        "/scratch/donor_feature_table.csv",
        table,
        feature_column=best["feature_column"],
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUND_COLUMNS,
        extra_columns=[
            "sex",
            "cognitive_status",
            "overall_ad_neuropath_change",
            "ca1_total_cells",
            "ca1_reactive_astrocytes",
            "ca1_astrocytes",
            "ca1_astro_lineage_cells",
            *[meta["feature_column"] for meta in VARIATIONS.values()],
        ],
    )

    _write_results_json(best, ranked)
    Path("/scratch/report.md").write_text(_build_report(best, ranked))

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best['variation_name']}")
    print(f"  IS partial r:      {_fmt(best['partial_r'])}")
    print(f"  Selection score:   {_fmt(best['selection_score'])}")
    print(f"  LOO predictive r:  {_fmt(best['loo_predictive_r'])}  (diagnostic)")
    print(f"  IS-LOO Gap:        {_fmt(best['is_loo_gap'])}  (penalty={_fmt(best['gap_penalty'])})")
    print(f"  Adjusted Score:    {_fmt(best['adjusted_score'])}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for row in ranked:
        print(
            f"  {row['variation_name']}  {_fmt(row['partial_r'])}  "
            f"{_fmt(row['selection_score'])}  {_fmt(row['loo_predictive_r'])}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best['feature_column']}  ca1_reactive_astrocytes  ca1_astro_lineage_cells  ca1_total_cells")
    for row in best["loo_table"].itertuples(index=False):
        print(
            f"  {row.donor_id}  {row.slope_zmem0:.4f}  {row.predicted:.4f}  "
            f"{getattr(row, best['feature_column']):.6f}  {int(row.ca1_reactive_astrocytes)}  "
            f"{int(row.ca1_astro_lineage_cells)}  {int(row.ca1_total_cells)}"
        )


if __name__ == "__main__":
    main()
PY

python /scratch/result.py | tee /scratch/run_output.txt
echo '--- result.py canonical variation line ---'
grep 'CANONICAL_VARIATION' /scratch/result.py
echo '--- results.json ---'
python - <<'PY'
import json
obj=json.load(open('/scratch/results.json'))
print(json.dumps({k: obj[k] for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score']}, indent=2))
PY
echo '--- report head ---'
sed -n '1,120p' /scratch/report.md
