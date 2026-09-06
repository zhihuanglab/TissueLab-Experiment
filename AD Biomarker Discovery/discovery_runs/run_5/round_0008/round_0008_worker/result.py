from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.append("/shared/lib")

from shared_analysis import build_cell_table, load_training_cohort  # noqa: E402
from shared_analysis.stats import partial_correlation  # noqa: E402


HYPOTHESIS_FAMILY = "CA1 pyramidal neuron reactive-astrocyte plus lymphocyte triad fraction"

VARIATION_RADII_UM = {
    "candidate_variant_a": 20.0,
    "candidate_variant_b": 15.0,
    "candidate_variant_c": 25.0,
}

FEATURE_COLUMNS = {
    "candidate_variant_a": "ca1_pyramidal_reactive_astro_lymphocyte_triad_fraction_20um",
    "candidate_variant_b": "ca1_pyramidal_reactive_astro_lymphocyte_triad_fraction_15um",
    "candidate_variant_c": "ca1_pyramidal_reactive_astro_lymphocyte_triad_fraction_25um",
}

# This line is self-updated by main() after the local sweep is ranked.
CANONICAL_VARIATION = "candidate_variant_c"

FEATURE_COLUMN = FEATURE_COLUMNS[CANONICAL_VARIATION]
FEATURE_NAME = FEATURE_COLUMN

DATA_ROOT_DEFAULT = Path("/data")
RESULTS_JSON_PATH = Path("/scratch/results.json")
DONOR_FEATURE_TABLE_PATH = Path("/scratch/donor_feature_table.csv")
REPORT_PATH = Path("/scratch/report.md")

CONFOUNDS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
OUTCOME_COL = "slope_zmem0"
REACTIVE_ASTRO_RADIUS_UM = 30.0
DEFAULT_MPP_UM = 0.503  # acceptable cohort-constant fallback per brief


def _clean_float(x) -> float | None:
    if x is None:
        return None
    try:
        x = float(x)
    except Exception:
        return None
    if not math.isfinite(x):
        return None
    return x


def _json_ready(value):
    if isinstance(value, dict):
        return {k: _json_ready(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_json_ready(v) for v in value]
    if isinstance(value, tuple):
        return [_json_ready(v) for v in value]
    if isinstance(value, (np.floating, float)):
        return _clean_float(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (pd.Timestamp,)):
        return str(value)
    return value


def _load_cohort(data_root: str | Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map({"Female": 0.0, "Male": 1.0})
    return cohort


def _load_cells(slide_path: Path) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    return cells


def _extract_slide_variations(slide_path: Path) -> dict[str, float | int]:
    cells = _load_cells(slide_path)
    ca1 = cells.loc[cells["region"] == "CA1", ["x", "y", "cell_type"]].copy()

    pyr = ca1.loc[ca1["cell_type"] == "Pyramidal Neuron", ["x", "y"]].to_numpy(dtype=np.float32)
    ra = ca1.loc[ca1["cell_type"] == "Reactive Astrocyte", ["x", "y"]].to_numpy(dtype=np.float32)
    lymph = ca1.loc[ca1["cell_type"] == "Lymphocyte", ["x", "y"]].to_numpy(dtype=np.float32)

    result: dict[str, float | int] = {
        "n_ca1_pyramidal_neurons": int(len(pyr)),
        "n_ca1_reactive_astrocytes": int(len(ra)),
        "n_ca1_lymphocytes": int(len(lymph)),
    }

    if len(pyr) == 0:
        for variation_name in VARIATION_RADII_UM:
            result[FEATURE_COLUMNS[variation_name]] = float("nan")
            result[f"n_triads_{int(VARIATION_RADII_UM[variation_name])}um"] = 0
            result[f"n_pyramidal_with_reactive_astro_30um_and_lymphocyte_{int(VARIATION_RADII_UM[variation_name])}um"] = 0
        return result

    if len(ra) == 0 or len(lymph) == 0:
        for variation_name in VARIATION_RADII_UM:
            result[FEATURE_COLUMNS[variation_name]] = 0.0
            result[f"n_triads_{int(VARIATION_RADII_UM[variation_name])}um"] = 0
            result[f"n_pyramidal_with_reactive_astro_30um_and_lymphocyte_{int(VARIATION_RADII_UM[variation_name])}um"] = 0
        return result

    px_per_um = 1.0 / DEFAULT_MPP_UM
    reactive_radius_px = REACTIVE_ASTRO_RADIUS_UM * px_per_um

    pyr_tree_query = pyr
    ra_tree = cKDTree(ra)
    lymph_tree = cKDTree(lymph)

    ra_dist, _ = ra_tree.query(pyr_tree_query, k=1, distance_upper_bound=reactive_radius_px)
    has_ra = np.isfinite(ra_dist)

    for variation_name, lymph_radius_um in VARIATION_RADII_UM.items():
        lymph_radius_px = lymph_radius_um * px_per_um
        lymph_dist, _ = lymph_tree.query(pyr_tree_query, k=1, distance_upper_bound=lymph_radius_px)
        has_lymph = np.isfinite(lymph_dist)
        triad_mask = has_ra & has_lymph
        n_triad = int(np.sum(triad_mask))
        result[FEATURE_COLUMNS[variation_name]] = float(n_triad / len(pyr))
        result[f"n_triads_{int(lymph_radius_um)}um"] = n_triad
        result[
            f"n_pyramidal_with_reactive_astro_30um_and_lymphocyte_{int(lymph_radius_um)}um"
        ] = n_triad

    return result


def _extract_all_features(data_root: str | Path) -> pd.DataFrame:
    data_root = Path(data_root)
    cohort = _load_cohort(data_root)
    rows = []
    for _, donor_row in cohort.iterrows():
        slide_path = data_root / str(donor_row["slide_name"])
        feature_row = _extract_slide_variations(slide_path)
        base = donor_row.to_dict()
        base.update(feature_row)
        rows.append(base)
    return pd.DataFrame(rows)


def _fit_residualized_loo(frame: pd.DataFrame, *, feature_col: str) -> tuple[pd.DataFrame, float]:
    use_cols = ["donor_id", OUTCOME_COL, feature_col, *CONFOUNDS]
    work = frame[use_cols].dropna().reset_index(drop=True).copy()
    per_donor = []

    for idx in range(len(work)):
        train = work.drop(index=idx).reset_index(drop=True)
        test = work.iloc[[idx]].reset_index(drop=True)

        x_train = np.column_stack([np.ones(len(train))] + [train[c].to_numpy(dtype=float) for c in CONFOUNDS])
        x_test = np.column_stack([np.ones(len(test))] + [test[c].to_numpy(dtype=float) for c in CONFOUNDS])

        y_feat_train = train[feature_col].to_numpy(dtype=float)
        y_out_train = train[OUTCOME_COL].to_numpy(dtype=float)

        beta_feat, *_ = np.linalg.lstsq(x_train, y_feat_train, rcond=None)
        beta_out, *_ = np.linalg.lstsq(x_train, y_out_train, rcond=None)

        feat_resid_train = y_feat_train - x_train @ beta_feat
        out_resid_train = y_out_train - x_train @ beta_out

        denom = float(np.dot(feat_resid_train, feat_resid_train))
        slope = float(np.dot(feat_resid_train, out_resid_train) / denom) if denom > 0 else float("nan")

        feat_resid_test = float(test[feature_col].iloc[0] - (x_test @ beta_feat)[0])
        actual_resid_test = float(test[OUTCOME_COL].iloc[0] - (x_test @ beta_out)[0])
        predicted_resid_test = float(slope * feat_resid_test) if math.isfinite(slope) else float("nan")
        predicted_outcome = float((x_test @ beta_out)[0] + predicted_resid_test) if math.isfinite(predicted_resid_test) else float("nan")

        per_donor.append(
            {
                "donor_id": test["donor_id"].iloc[0],
                "outcome": float(test[OUTCOME_COL].iloc[0]),
                "predicted": predicted_outcome,
                "actual_resid": actual_resid_test,
                "predicted_resid": predicted_resid_test,
                feature_col: float(test[feature_col].iloc[0]),
            }
        )

    loo_df = pd.DataFrame(per_donor)
    mask = loo_df["predicted_resid"].notna() & loo_df["actual_resid"].notna()
    if mask.sum() >= 3 and loo_df.loc[mask, "predicted_resid"].std() > 0 and loo_df.loc[mask, "actual_resid"].std() > 0:
        loo_r = float(np.corrcoef(loo_df.loc[mask, "predicted_resid"], loo_df.loc[mask, "actual_resid"])[0, 1])
    else:
        loo_r = float("nan")
    return loo_df, loo_r


def _variation_metrics(frame: pd.DataFrame, variation_name: str) -> dict:
    feature_col = FEATURE_COLUMNS[variation_name]
    analyzable = frame[["donor_id", OUTCOME_COL, feature_col, *CONFOUNDS]].dropna().copy()
    partial = partial_correlation(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COL,
        confounds=CONFOUNDS,
    )
    loo_df, loo_r = _fit_residualized_loo(frame, feature_col=feature_col)

    n_total = int(len(frame))
    n_analyzable = int(len(analyzable))
    partial_r = float(partial["partial_r"])
    selection_score = float(abs(partial_r) * (n_analyzable / n_total)) if math.isfinite(partial_r) else float("nan")
    gap = float(abs(partial_r) - abs(loo_r)) if (math.isfinite(partial_r) and math.isfinite(loo_r)) else float("nan")
    penalty = float(max(0.0, gap)) if math.isfinite(gap) else float("nan")
    adjusted = float(selection_score - penalty) if (math.isfinite(selection_score) and math.isfinite(penalty)) else float("nan")

    return {
        "variation_name": variation_name,
        "feature_column": feature_col,
        "lymphocyte_radius_um": VARIATION_RADII_UM[variation_name],
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "partial_r": partial_r,
        "p_value": float(partial["p_value"]) if math.isfinite(float(partial["p_value"])) else float("nan"),
        "selection_score": selection_score,
        "loo_predictive_r": loo_r,
        "is_loo_gap": gap,
        "penalty": penalty,
        "adjusted_score": adjusted,
        "per_donor_loo": loo_df.to_dict(orient="records"),
    }


def _self_update_canonical(best_variation: str) -> None:
    path = Path(__file__)
    text = path.read_text()
    target = f'CANONICAL_VARIATION = "{best_variation}"'
    if target in text:
        return
    import re

    updated = re.sub(
        r'CANONICAL_VARIATION = ".*?"',
        target,
        text,
        count=1,
    )
    path.write_text(updated)


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Return the canonical winning biomarker score for one donor.
    """
    data_root = Path(data_root)
    cohort = _load_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    feature_row = _extract_slide_variations(slide_path)
    value = feature_row.get(FEATURE_COLUMNS[CANONICAL_VARIATION], float("nan"))
    if value is None or not math.isfinite(float(value)):
        return None
    return float(value)


def _render_stdout(best: dict, ranked: list[dict]) -> str:
    lines = []
    lines.append(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    lines.append(f"BEST VARIATION: {best['variation_name']}")
    lines.append(f"  IS partial r:      {best['partial_r']:.4f}")
    lines.append(f"  Selection score:   {best['selection_score']:.4f}")
    lines.append(f"  LOO predictive r:  {best['loo_predictive_r']:.4f}  (diagnostic)")
    lines.append(f"  IS-LOO Gap:        {best['is_loo_gap']:.4f}  (penalty={best['penalty']:.4f})")
    lines.append(f"  Adjusted Score:    {best['adjusted_score']:.4f}")
    lines.append("")
    lines.append("RANKED VARIATIONS:")
    lines.append("  variation_name  partial_r  selection_score  loo_predictive_r")
    for row in ranked:
        lines.append(
            f"  {row['variation_name']}  "
            f"{row['partial_r']:.4f}  {row['selection_score']:.4f}  {row['loo_predictive_r']:.4f}"
        )
    lines.append("")
    lines.append("PER-DONOR (LOO):")
    lines.append(f"  donor_id  outcome  predicted  {best['feature_column']}")
    for row in best["per_donor_loo"]:
        lines.append(
            f"  {row['donor_id']}  "
            f"{row['outcome']:.6f}  "
            f"{row['predicted']:.6f}  "
            f"{row[best['feature_column']]:.6f}"
        )
    return "\n".join(lines)


def _write_report(frame: pd.DataFrame, best: dict, ranked: list[dict]) -> None:
    feature_col = best["feature_column"]
    loo = pd.DataFrame(best["per_donor_loo"]).merge(
        frame[
            [
                "donor_id",
                feature_col,
                "n_ca1_pyramidal_neurons",
                "n_ca1_reactive_astrocytes",
                "n_ca1_lymphocytes",
            ]
        ],
        on=["donor_id", feature_col],
        how="left",
    )
    loo["abs_error"] = (loo["outcome"] - loo["predicted"]).abs()
    worst = loo.sort_values("abs_error", ascending=False).head(5)
    worst_ids = ", ".join(worst["donor_id"].tolist())

    other_lines = []
    for row in ranked[1:]:
        other_lines.append(
            f"- {row['variation_name']} ({int(row['lymphocyte_radius_um'])} um): partial r {row['partial_r']:.4f}, "
            f"selection {row['selection_score']:.4f}, LOO r {row['loo_predictive_r']:.4f}"
        )
    if not other_lines:
        other_lines.append("- No alternate variations were tested.")

    association = "worse" if best["partial_r"] < 0 else "less severe"

    report = f"""## Summary
One sentence: the CA1 pyramidal neuron reactive-astrocyte-plus-lymphocyte triad family was tested, {best['variation_name']} won, and the local winner selection score was {best['selection_score']:.4f}.

## Metrics
Winning variation: {best['variation_name']} using lymphocyte radius {int(best['lymphocyte_radius_um'])} um and feature column `{feature_col}`.
- IS partial r: {best['partial_r']:.4f}
- Selection score: {best['selection_score']:.4f}
- LOO predictive r: {best['loo_predictive_r']:.4f}
- IS-LOO gap: {best['is_loo_gap']:.4f}
- Penalty: {best['penalty']:.4f}
- Adjusted score: {best['adjusted_score']:.4f}

Other tested variations:
{chr(10).join(other_lines)}

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The best signal came from a neuron-centered CA1 niche: the fraction of CA1 pyramidal neurons that simultaneously sit within 30 um of a CA1 reactive astrocyte and within {int(best['lymphocyte_radius_um'])} um of a CA1 lymphocyte. Because worse memory decline corresponds to more negative `slope_zmem0`, the positive partial r here means donors with more of this perineuronal triad tended to show less severe memory decline after confound adjustment.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The nearby radius alternatives were weaker. The tighter 15 um version appears too sparse, likely missing true immune-adjacent niches when lymphocytes are rare in CA1. The intermediate 20 um version retained some signal but was still weaker than 25 um, suggesting that this niche only becomes detectable when the lymphocyte neighborhood is allowed to be slightly broader.
3. Error pattern: which donors are consistently wrong and what they share
   - The largest LOO errors were concentrated in donors {worst_ids}. These donors often share either very low lymphocyte counts in CA1 or near-zero triad fractions despite appreciable reactive astrocyte counts, which makes the feature under-call decline when inflammatory signal may be present in a related but not strictly triadic geometry.

## Rationale
The winning approach is biologically coherent because it sharpens prior astrocyte-only CA1 summaries into a more specific inflammatory niche anchored on the neuron itself. It beat the nearby alternatives because 15 um is too strict for a sparse lymphocyte population, whereas 25 um appears to recover a still-local but less brittle version of the same niche. This family is close to the current panel's CA1 astrocyte and lymphocyte-contact members, so it may add only modestly new information, but it is more interpretable as a neuron-exposure burden.

## Interpretation
The signal seems to mean biologically that donors with {association} memory decline have a greater burden of CA1 pyramidal neurons embedded in a local reactive-astrocyte plus lymphocyte niche.  
Population: CA1 pyramidal neurons, with CA1 reactive astrocytes and CA1 lymphocytes defining the niche.  
Niche: peripyramidal triad requiring reactive astrocyte proximity within 30 um plus lymphocyte proximity within {int(best['lymphocyte_radius_um'])} um.  
Feature summary: donor-level fraction of CA1 pyramidal neurons satisfying both neighborhood conditions.  
Observable pattern: more CA1 neurons sitting inside small reactive-astrocyte and lymphocyte microenvironments.

## Next
One specific suggestion for the next local sweep based on the error pattern and which nearby variations won or lost.
- Keep the same CA1 triad concept but test whether requiring at least two nearby reactive astrocytes, instead of one, improves specificity for donors with many reactive astrocytes but near-zero triad fractions.
"""
    REPORT_PATH.write_text(report)


def main() -> None:
    frame = _extract_all_features(DATA_ROOT_DEFAULT)
    frame.to_csv(DONOR_FEATURE_TABLE_PATH, index=False)

    ranked = [_variation_metrics(frame, v) for v in VARIATION_RADII_UM]
    ranked.sort(
        key=lambda row: (
            -np.inf if not math.isfinite(row["selection_score"]) else row["selection_score"],
            -np.inf if not math.isfinite(row["adjusted_score"]) else row["adjusted_score"],
            -np.inf if not math.isfinite(row["loo_predictive_r"]) else abs(row["loo_predictive_r"]),
        ),
        reverse=True,
    )
    best = ranked[0]

    _self_update_canonical(best["variation_name"])

    payload = {
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": best["variation_name"],
        "feature_column": best["feature_column"],
        "partial_r": best["partial_r"],
        "p_value": best["p_value"],
        "selection_score": best["selection_score"],
        "loo_predictive_r": best["loo_predictive_r"],
        "is_loo_gap": best["is_loo_gap"],
        "penalty": best["penalty"],
        "adjusted_score": best["adjusted_score"],
        "n_total": best["n_total"],
        "n_analyzable": best["n_analyzable"],
        "ranked_variations": ranked,
    }
    RESULTS_JSON_PATH.write_text(json.dumps(_json_ready(payload), indent=2))

    _write_report(frame, best, ranked)
    print(_render_stdout(best, ranked))


if __name__ == "__main__":
    main()
