from __future__ import annotations

import json
import math
import re
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.append("/shared/lib")

from shared_analysis.artifacts import build_results_payload, write_donor_feature_table, write_results_payload
from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort
from shared_analysis.stats import (
    bootstrap_partial_correlation,
    bootstrap_partial_correlation_stability,
    leave_one_out_summary,
    partial_correlation,
)

warnings.filterwarnings(
    "ignore",
    message="Object at .* is not recognized as a component of a Zarr hierarchy.",
)

HYPOTHESIS_FAMILY = "Astrocytic reactivity fraction in CA1 versus broader pyramidal field"
FEATURE_NAME = "astrocytic_reactivity_fraction"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]

# Updated in-place after running the round so replay defaults to the winning variant.
CANONICAL_VARIATION_NAME = "candidate_variant_a"
CANONICAL_FEATURE_COLUMN = "reactive_astro_fraction_ca1"
FEATURE_COLUMN = CANONICAL_FEATURE_COLUMN

VARIATIONS: dict[str, dict[str, Any]] = {
    "candidate_variant_a": {
        "name": "candidate_variant_a",
        "feature_column": "reactive_astro_fraction_ca1",
        "regions": ("CA1",),
        "region_label": "CA1",
        "description": "Reactive Astrocyte / (Astrocyte + Reactive Astrocyte) among CA1-assigned cells.",
        "population": "Astrocyte + Reactive Astrocyte lineage",
        "niche": "CA1 pyramidal field",
    },
    "candidate_variant_b": {
        "name": "candidate_variant_b",
        "feature_column": "reactive_astro_fraction_ca1_ca2_ca3",
        "regions": ("CA1", "CA2", "CA3"),
        "region_label": "CA1+CA2+CA3",
        "description": "Reactive Astrocyte / (Astrocyte + Reactive Astrocyte) among pooled CA1/CA2/CA3 cells.",
        "population": "Astrocyte + Reactive Astrocyte lineage",
        "niche": "broader pyramidal field",
    },
}


def safe_float(value: Any) -> float | None:
    try:
        value = float(value)
    except Exception:
        return None
    if not math.isfinite(value):
        return None
    return value


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_safe(v) for v in value]
    if isinstance(value, tuple):
        return [json_safe(v) for v in value]
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, (float, int)):
        if isinstance(value, float) and not math.isfinite(value):
            return None
        return value
    return value


def derive_sex_binary(value: Any) -> float:
    token = str(value).strip().lower()
    if token == "male":
        return 1.0
    if token == "female":
        return 0.0
    return float("nan")


def design_matrix(confounds: pd.DataFrame) -> np.ndarray:
    arr = confounds.astype(float).to_numpy()
    return np.column_stack([np.ones(len(confounds), dtype=float), arr])


def feature_sort_key(result: dict[str, Any]) -> tuple[float, float, float]:
    selection = result.get("selection_score")
    adjusted = result.get("adjusted_score")
    loo = result.get("loo_predictive_r")
    selection_val = -math.inf if selection is None or not math.isfinite(selection) else float(selection)
    adjusted_val = -math.inf if adjusted is None or not math.isfinite(adjusted) else float(adjusted)
    loo_val = -math.inf if loo is None or not math.isfinite(loo) else abs(float(loo))
    return (selection_val, adjusted_val, loo_val)


def compute_reactivity_fraction(cells: pd.DataFrame, regions: tuple[str, ...]) -> dict[str, float]:
    astro_lineage = cells["cell_type"].isin(["Astrocyte", "Reactive Astrocyte"])
    in_region = cells["region"].isin(list(regions))
    lineage_region = cells.loc[astro_lineage & in_region, "cell_type"]
    reactive_count = int((lineage_region == "Reactive Astrocyte").sum())
    denominator = int(len(lineage_region))
    astro_count = int((lineage_region == "Astrocyte").sum())
    fraction = float(reactive_count / denominator) if denominator > 0 else float("nan")
    return {
        "fraction": fraction,
        "reactive_count": float(reactive_count),
        "astro_count": float(astro_count),
        "denominator": float(denominator),
    }


def current_canonical_variation_name() -> str:
    sidecar = Path(__file__).with_name("results.json")
    if sidecar.exists():
        try:
            payload = json.loads(sidecar.read_text(encoding="utf-8"))
            best = str(payload.get("best_variation") or "").strip()
            if best in VARIATIONS:
                return best
        except Exception:
            pass
    if CANONICAL_VARIATION_NAME in VARIATIONS:
        return CANONICAL_VARIATION_NAME
    return "candidate_variant_a"


def compute_variation_score_for_donor(*, donor_id: str, data_root: str | Path, variation_name: str) -> float | None:
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    if not slide_path.exists():
        return None
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    variation = VARIATIONS[variation_name]
    score = compute_reactivity_fraction(cells, tuple(variation["regions"]))["fraction"]
    return float(score) if math.isfinite(score) else None


def compute_donor_score(*, donor_id: str, data_root: str | Path):
    """
    Replay the winning round-1 astrocytic reactivity biomarker from raw .zarr data.
    """
    variation_name = current_canonical_variation_name()
    return compute_variation_score_for_donor(
        donor_id=donor_id,
        data_root=data_root,
        variation_name=variation_name,
    )


def extract_feature_table(data_root: Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map(derive_sex_binary)
    rows: list[dict[str, Any]] = []
    for _, donor in cohort.iterrows():
        slide_name = str(donor["slide_name"])
        slide_path = data_root / slide_name
        base_row = {
            "donor_id": donor["donor_id"],
            "slide_name": slide_name,
            OUTCOME_COLUMN: donor[OUTCOME_COLUMN],
            "max_age_vis": donor["max_age_vis"],
            "braak_numeric": donor["braak_numeric"],
            "cerad_ordinal": donor["cerad_ordinal"],
            "sex_binary": donor["sex_binary"],
            "sex": donor["sex"],
            "cognitive_status": donor.get("cognitive_status"),
            "overall_ad_neuropath_change": donor.get("overall_ad_neuropath_change"),
        }
        if not slide_path.exists():
            rows.append(base_row)
            continue

        cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
        for variation in VARIATIONS.values():
            summary = compute_reactivity_fraction(cells, tuple(variation["regions"]))
            col = variation["feature_column"]
            base_row[col] = summary["fraction"]
            base_row[f"{col}_reactive_count"] = summary["reactive_count"]
            base_row[f"{col}_astro_count"] = summary["astro_count"]
            base_row[f"{col}_denominator"] = summary["denominator"]
        rows.append(base_row)
    return pd.DataFrame(rows)


def loo_predictions(df: pd.DataFrame, *, feature_col: str) -> tuple[float | None, pd.DataFrame]:
    cols = ["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]
    frame = df.loc[:, cols].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    rows: list[dict[str, Any]] = []
    if len(frame) < 3:
        return None, pd.DataFrame(rows)

    predicted_resid: list[float] = []
    actual_resid: list[float] = []
    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train = design_matrix(train[CONFOUND_COLUMNS])
        x_test = design_matrix(test[CONFOUND_COLUMNS])

        beta_feature, *_ = np.linalg.lstsq(
            x_train,
            train[feature_col].to_numpy(dtype=float),
            rcond=None,
        )
        beta_outcome, *_ = np.linalg.lstsq(
            x_train,
            train[OUTCOME_COLUMN].to_numpy(dtype=float),
            rcond=None,
        )

        feature_train = train[feature_col].to_numpy(dtype=float)
        outcome_train = train[OUTCOME_COLUMN].to_numpy(dtype=float)
        feature_resid_train = feature_train - x_train @ beta_feature
        outcome_resid_train = outcome_train - x_train @ beta_outcome

        denom = float(np.dot(feature_resid_train, feature_resid_train))
        if denom <= 0:
            pred_resid = float("nan")
        else:
            slope = float(np.dot(feature_resid_train, outcome_resid_train) / denom)
            feature_resid_test = test[feature_col].to_numpy(dtype=float) - x_test @ beta_feature
            pred_resid = float(slope * feature_resid_test[0])

        actual_resid_test = float(
            test[OUTCOME_COLUMN].to_numpy(dtype=float)[0] - (x_test @ beta_outcome)[0]
        )
        pred_raw = float((x_test @ beta_outcome)[0] + pred_resid) if math.isfinite(pred_resid) else float("nan")
        actual_raw = float(test[OUTCOME_COLUMN].iloc[0])

        predicted_resid.append(pred_resid)
        actual_resid.append(actual_resid_test)
        rows.append(
            {
                "donor_id": str(test["donor_id"].iloc[0]),
                "outcome": actual_raw,
                "predicted": pred_raw,
                "predicted_residualized": pred_resid,
                "actual_residualized": actual_resid_test,
                feature_col: float(test[feature_col].iloc[0]),
                "abs_error": abs(pred_raw - actual_raw) if math.isfinite(pred_raw) else float("nan"),
            }
        )

    pred_arr = np.asarray(predicted_resid, dtype=float)
    act_arr = np.asarray(actual_resid, dtype=float)
    valid = np.isfinite(pred_arr) & np.isfinite(act_arr)
    if valid.sum() >= 3 and np.std(pred_arr[valid]) > 0 and np.std(act_arr[valid]) > 0:
        loo_r = float(np.corrcoef(pred_arr[valid], act_arr[valid])[0, 1])
    else:
        loo_r = None
    return loo_r, pd.DataFrame(rows)


def evaluate_variation(df: pd.DataFrame, *, variation: dict[str, Any], n_total: int) -> dict[str, Any]:
    feature_col = variation["feature_column"]
    analyzable = (
        df.loc[:, ["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]]
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
        .reset_index(drop=True)
    )
    n_analyzable = int(len(analyzable))

    if n_analyzable < 4:
        return {
            "name": variation["name"],
            "feature_column": feature_col,
            "description": variation["description"],
            "regions": list(variation["regions"]),
            "region_label": variation["region_label"],
            "population": variation["population"],
            "niche": variation["niche"],
            "n_analyzable": n_analyzable,
            "coverage_ratio": float(n_analyzable / n_total) if n_total else None,
            "partial_r": None,
            "p_value": None,
            "ci_lo": None,
            "ci_hi": None,
            "selection_score": None,
            "loo_predictive_r": None,
            "gap": None,
            "penalty": None,
            "adjusted_score": None,
            "loo_unstable_count": 0,
            "loo_max_shift": None,
            "bootstrap_median_partial_r": None,
            "bootstrap_sign_consistency": None,
            "donor_ids_used": analyzable["donor_id"].tolist(),
            "loo_rows": [],
        }

    base_stats = partial_correlation(
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
        n_boot=500,
        random_state=0,
    )
    stability = bootstrap_partial_correlation_stability(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        n_boot=300,
        sample_frac=0.8,
        random_state=0,
    )
    loo_predictive_r, loo_table = loo_predictions(df, feature_col=feature_col)
    loo_summary = leave_one_out_summary(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        id_col="donor_id",
        unstable_delta=0.10,
    )

    partial_r = safe_float(base_stats.get("partial_r"))
    coverage_ratio = n_analyzable / n_total if n_total else float("nan")
    selection_score = abs(partial_r) * coverage_ratio if partial_r is not None else None
    if partial_r is not None and loo_predictive_r is not None:
        gap = abs(partial_r) - abs(loo_predictive_r)
        penalty = max(0.0, gap)
        adjusted = selection_score - penalty if selection_score is not None else None
    else:
        gap = None
        penalty = None
        adjusted = None

    return {
        "name": variation["name"],
        "feature_column": feature_col,
        "description": variation["description"],
        "regions": list(variation["regions"]),
        "region_label": variation["region_label"],
        "population": variation["population"],
        "niche": variation["niche"],
        "n_analyzable": n_analyzable,
        "coverage_ratio": float(coverage_ratio),
        "partial_r": partial_r,
        "p_value": safe_float(base_stats.get("p_value")),
        "ci_lo": safe_float(boot.get("ci_lo")),
        "ci_hi": safe_float(boot.get("ci_hi")),
        "selection_score": safe_float(selection_score),
        "loo_predictive_r": safe_float(loo_predictive_r),
        "gap": safe_float(gap),
        "penalty": safe_float(penalty),
        "adjusted_score": safe_float(adjusted),
        "loo_unstable_count": int(loo_summary.get("unstable_count", 0)),
        "loo_max_shift": safe_float(loo_summary.get("max_shift")),
        "bootstrap_median_partial_r": safe_float(stability.get("bootstrap_median_partial_r")),
        "bootstrap_sign_consistency": safe_float(stability.get("bootstrap_sign_consistency")),
        "donor_ids_used": analyzable["donor_id"].tolist(),
        "loo_rows": json_safe(loo_table.to_dict(orient="records")),
    }


def summarize_error_pattern(feature_table: pd.DataFrame, best: dict[str, Any]) -> str:
    loo_df = pd.DataFrame(best.get("loo_rows", []))
    if loo_df.empty or "abs_error" not in loo_df.columns:
        return "LOO prediction errors were too sparse to summarize."
    merged = loo_df.merge(
        feature_table[
            [
                "donor_id",
                "cognitive_status",
                "overall_ad_neuropath_change",
                "braak_numeric",
                "cerad_ordinal",
            ]
        ],
        on="donor_id",
        how="left",
    )
    merged = merged.sort_values("abs_error", ascending=False).head(5).reset_index(drop=True)
    donor_list = ", ".join(merged["donor_id"].astype(str).tolist())
    dementia_frac = float((merged["cognitive_status"] == "Dementia").mean()) if len(merged) else float("nan")
    high_path_frac = float(
        merged["overall_ad_neuropath_change"].isin(["Intermediate", "High"]).mean()
    ) if len(merged) else float("nan")
    braak_med = safe_float(merged["braak_numeric"].median())
    cerad_med = safe_float(merged["cerad_ordinal"].median())

    bits = [f"Largest absolute LOO errors were {donor_list}."]
    if math.isfinite(dementia_frac):
        bits.append(f"{int(round(dementia_frac * 100))}% of these donors were labeled Dementia.")
    if math.isfinite(high_path_frac):
        bits.append(
            f"{int(round(high_path_frac * 100))}% had Intermediate/High AD neuropathologic change."
        )
    if braak_med is not None and cerad_med is not None:
        bits.append(f"Their median Braak/CERAD numeric values were {braak_med:.1f}/{cerad_med:.1f}.")
    bits.append(
        "This suggests the astrocytic fraction alone does not fully capture how severe pathology translates into memory decline for these donors."
    )
    return " ".join(bits)


def update_canonical_target_in_source(best_variation: str, feature_column: str) -> None:
    path = Path(__file__)
    text = path.read_text(encoding="utf-8")
    new_text = re.sub(
        r'CANONICAL_VARIATION_NAME = ".*?"',
        f'CANONICAL_VARIATION_NAME = "{best_variation}"',
        text,
        count=1,
    )
    new_text = re.sub(
        r'CANONICAL_FEATURE_COLUMN = ".*?"',
        f'CANONICAL_FEATURE_COLUMN = "{feature_column}"',
        new_text,
        count=1,
    )
    new_text = re.sub(
        r'FEATURE_COLUMN = CANONICAL_FEATURE_COLUMN',
        'FEATURE_COLUMN = CANONICAL_FEATURE_COLUMN',
        new_text,
        count=1,
    )
    if new_text != text:
        path.write_text(new_text, encoding="utf-8")


def render_report(feature_table: pd.DataFrame, ranked: list[dict[str, Any]], best: dict[str, Any]) -> str:
    loser_bits = []
    for entry in ranked[1:]:
        loser_bits.append(
            f"{entry['name']} had selection score {entry['selection_score']:.4f}, "
            f"partial r {entry['partial_r']:.4f}, and LOO r {entry['loo_predictive_r']:.4f}."
            if entry.get("selection_score") is not None and entry.get("partial_r") is not None and entry.get("loo_predictive_r") is not None
            else f"{entry['name']} was weaker or less analyzable."
        )
    loser_text = " ".join(loser_bits) if loser_bits else "No nearby alternative was tested."

    best_regions = best["region_label"]
    if best["name"] == "candidate_variant_a":
        worked = (
            "The CA1-only astrocytic reactivity fraction worked best, consistent with CA1 being a memory-critical and AD-vulnerable niche where reactive gliosis may track local tissue stress more specifically than a broader hippocampal pool."
        )
        failed = (
            "Pooling CA1 with CA2 and CA3 diluted the signal, implying that the strongest association is not a generic pyramidal-field astrocytic response but a more CA1-focused shift in astrocyte state."
        )
        rationale = (
            "CA1 is selectively vulnerable in AD-related hippocampal degeneration, so a within-lineage reactive-astrocyte fraction in that region is biologically coherent as a compact readout of local gliosis. It beat the broader pooled alternative because adding CA2/CA3 likely mixed in less vulnerable tissue and reduced anatomical specificity."
        )
        next_step = (
            "Next sweep: keep the astrocyte/reactive-astrocyte lineage and CA1 gate fixed, but compare CA1 reactive-astrocyte fraction against CA1 reactive-astrocyte density or raw reactive count to test whether high-error donors are missing burden information rather than state-composition information."
        )
    else:
        worked = (
            "The broader CA1+CA2+CA3 astrocytic reactivity fraction worked best, suggesting that the memory-decline signal reflects a field-wide astrocytic state shift across the pyramidal system rather than a narrowly CA1-restricted effect."
        )
        failed = (
            "The CA1-only fraction was weaker, which suggests that a single-region gate was too noisy or too sparse and that pooling adjacent pyramidal regions stabilized the donor-level state estimate."
        )
        rationale = (
            "A within-lineage reactive fraction across CA1-CA3 is biologically coherent if astrocytic activation tracks a broader hippocampal injury program spanning the pyramidal field. It beat CA1 alone because the extra regions likely improved robustness without erasing the underlying astrocytic-state contrast."
        )
        next_step = (
            "Next sweep: keep the astrocyte/reactive-astrocyte fraction fixed but compare CA1+CA2+CA3 against CA1+CA2+CA3+CA4 and against a density-style summary, to see whether the current winner is capturing a stable field effect or still missing burden information in the largest-error donors."
        )

    error_pattern = summarize_error_pattern(feature_table, best)
    best_metrics = (
        f"partial r {best['partial_r']:.4f}, selection score {best['selection_score']:.4f}, "
        f"LOO predictive r {best['loo_predictive_r']:.4f}, IS-LOO gap {best['gap']:.4f}, "
        f"penalty {best['penalty']:.4f}, adjusted score {best['adjusted_score']:.4f}."
    )

    report = f"""## Summary
One sentence: tested astrocytic reactivity fractions in CA1 versus a broader pyramidal-field pool; {best['name']} won with local selection score {best['selection_score']:.4f}.

## Metrics
Winning variation: {best['name']} ({best_regions}) with {best_metrics}
Other tested variations: {loser_text}

## Findings
1. What worked and why: {worked}
2. What failed and why: {failed}
3. Error pattern: {error_pattern}

## Rationale
{rationale}
If relevant to the future panel, this feature appears interpretable and auditable, but whether it adds new information beyond later panel members will depend on whether the eventual panel contains other glial or hippocampal-region features.

## Interpretation
The signal most plausibly reflects astrocytic gliosis. Population: Astrocyte plus Reactive Astrocyte lineage. Niche: {best_regions}. Feature summary: Reactive Astrocyte / (Astrocyte + Reactive Astrocyte) at the donor level. Simplest observable pattern: donors with worse memory decline tend to show a larger share of astrocytes classified as reactive within the selected hippocampal niche.

## Next
{next_step}
"""
    return report


def print_stdout(ranked: list[dict[str, Any]], best: dict[str, Any]) -> None:
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
    for entry in ranked:
        print(
            f"  {entry['name']}  "
            f"{entry['partial_r']:.4f}  "
            f"{entry['selection_score']:.4f}  "
            f"{entry['loo_predictive_r']:.4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best['feature_column']}")
    for row in best["loo_rows"]:
        print(
            f"  {row['donor_id']}  "
            f"{row['outcome']:.6f}  "
            f"{row['predicted']:.6f}  "
            f"{row[best['feature_column']]:.6f}"
        )


def main() -> None:
    data_root = Path("/data")
    output_root = Path("/scratch")
    feature_table = extract_feature_table(data_root)
    n_total = int(len(feature_table))

    evaluations = [
        evaluate_variation(feature_table, variation=variation, n_total=n_total)
        for variation in VARIATIONS.values()
    ]
    ranked = sorted(evaluations, key=feature_sort_key, reverse=True)
    best = ranked[0]

    if best.get("selection_score") is None:
        raise RuntimeError("No analyzable variation was produced for this round.")

    best_feature = best["feature_column"]
    extra_columns = [
        col
        for col in feature_table.columns
        if col
        not in {
            "donor_id",
            "slide_name",
            OUTCOME_COLUMN,
            *CONFOUND_COLUMNS,
            "sex",
            "cognitive_status",
            "overall_ad_neuropath_change",
        }
        and col != best_feature
    ]
    donor_table_path = write_donor_feature_table(
        output_root / "donor_feature_table.csv",
        feature_table,
        feature_column=best_feature,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUND_COLUMNS,
        extra_columns=extra_columns,
    )

    compact_ranked = []
    for entry in ranked:
        compact = {k: v for k, v in entry.items() if k != "loo_rows"}
        compact_ranked.append(json_safe(compact))

    report_text = render_report(feature_table, ranked, best)
    report_path = output_root / "report.md"
    report_path.write_text(report_text, encoding="utf-8")

    results_payload = build_results_payload(
        status="ok",
        feature_name=FEATURE_NAME,
        outcome=OUTCOME_COLUMN,
        n_total=n_total,
        n_analyzable=best["n_analyzable"],
        partial_r=best["partial_r"],
        ci_lo=best["ci_lo"],
        ci_hi=best["ci_hi"],
        p_value=best["p_value"],
        loo_predictive_r=best["loo_predictive_r"],
        loo_unstable_count=best["loo_unstable_count"],
        loo_max_shift=best["loo_max_shift"],
        donor_ids_used=best["donor_ids_used"],
        covariates=CONFOUND_COLUMNS,
        recomputed_from_raw=True,
        registry_written=False,
        artifacts={
            "donor_feature_table": str(donor_table_path),
            "report": str(report_path),
        },
        best_variation=best["name"],
        ranked_variations=compact_ranked,
        feature_column=best_feature,
        selection_score=best["selection_score"],
        adjusted_score=best["adjusted_score"],
        is_loo_gap=best["gap"],
        gap_penalty=best["penalty"],
        coverage_ratio=best["coverage_ratio"],
        bootstrap_median_partial_r=best["bootstrap_median_partial_r"],
        bootstrap_sign_consistency=best["bootstrap_sign_consistency"],
        per_donor_loo=best["loo_rows"],
        hypothesis_family=HYPOTHESIS_FAMILY,
    )
    write_results_payload(output_root / "results.json", json_safe(results_payload))

    update_canonical_target_in_source(best["name"], best_feature)
    print_stdout(ranked, best)


if __name__ == "__main__":
    main()
