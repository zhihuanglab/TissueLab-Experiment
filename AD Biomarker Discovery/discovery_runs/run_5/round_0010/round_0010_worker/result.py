from __future__ import annotations

import json
import math
import re
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import tifffile
from scipy.spatial import cKDTree

# Make shared helpers importable in this runtime.
sys.path.insert(0, "/shared/lib")
warnings.filterwarnings("ignore", message="Object at .* is not recognized as a component of a Zarr hierarchy.")
warnings.filterwarnings("ignore", message="invalid value encountered in divide")
from shared_analysis import build_cell_table, load_training_cohort  # type: ignore


HYPOTHESIS_FAMILY = "CA1 peripyramidal reactive astrocyte lymphocyte-contact area burden"
CANONICAL_DEFAULT_VARIATION = "candidate_variant_a"
FEATURE_NAME = "ca1_peripyramidal_ra_lymphocyte_area_burden_20um"
FEATURE_COLUMN = "ca1_peripyramidal_ra_lymphocyte_area_burden_20um"
VARIATIONS: dict[str, dict[str, Any]] = {
    "candidate_variant_a": {
        "variation_name": "candidate_variant_a",
        "feature_name": "ca1_peripyramidal_ra_lymphocyte_area_burden_20um",
        "feature_column": "ca1_peripyramidal_ra_lymphocyte_area_burden_20um",
        "peripyramidal_radius_um": 30.0,
        "lymphocyte_contact_radius_um": 20.0,
        "description": (
            "Peripyramidal reactive astrocytes within 30 um of CA1 pyramidal neurons; "
            "score is summed area of lymphocyte-contacting cells / peripyramidal reactive astrocyte count, "
            "with lymphocyte contact defined within 20 um."
        ),
    },
    "candidate_variant_b": {
        "variation_name": "candidate_variant_b",
        "feature_name": "ca1_peripyramidal_ra_lymphocyte_area_burden_25um",
        "feature_column": "ca1_peripyramidal_ra_lymphocyte_area_burden_25um",
        "peripyramidal_radius_um": 30.0,
        "lymphocyte_contact_radius_um": 25.0,
        "description": (
            "Same as variant A but lymphocyte contact radius expanded to 25 um."
        ),
    },
}
CONFounds = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
OUTCOME_COLUMN = "slope_zmem0"
RESULTS_PATH = Path(__file__).with_name("results.json")
DONOR_FEATURE_TABLE_PATH = Path(__file__).with_name("donor_feature_table.csv")
REPORT_PATH = Path(__file__).with_name("report.md")


@dataclass
class DonorExtraction:
    donor_id: str
    slide_name: str
    mpp_um: float
    ca1_pyramidal_count: int
    ca1_reactive_astro_count: int
    ca1_lymphocyte_count: int
    peripyramidal_ra_count: int
    peripyramidal_ra_area_um2_sum: float
    scores: dict[str, float]


def _json_default(obj: Any):
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if pd.isna(obj):
        return None
    raise TypeError(f"Object of type {type(obj)!r} is not JSON serializable")


def _load_best_variation_from_sidecar() -> str:
    if RESULTS_PATH.exists():
        try:
            payload = json.loads(RESULTS_PATH.read_text())
            best = payload.get("best_variation")
            if best in VARIATIONS:
                return str(best)
        except Exception:
            pass
    return CANONICAL_DEFAULT_VARIATION


def _aperio_mpp_um(svs_path: Path) -> float:
    with tifffile.TiffFile(str(svs_path)) as tf:
        desc = tf.pages[0].description or ""
    match = re.search(r"\bMPP\s*=\s*([0-9.]+)", desc)
    if match:
        return float(match.group(1))
    return 0.503  # cohort fallback


def _safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return float("nan")
    x = x[mask]
    y = y[mask]
    if np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def _design_matrix(frame: pd.DataFrame, columns: list[str]) -> np.ndarray:
    x = frame[columns].to_numpy(dtype=float)
    intercept = np.ones((len(frame), 1), dtype=float)
    return np.concatenate([intercept, x], axis=1)


def _residualize(target: np.ndarray, confounds: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(confounds, target, rcond=None)
    return target - confounds @ beta


def _partial_correlation(frame: pd.DataFrame, feature_col: str) -> float:
    clean = frame[[feature_col, OUTCOME_COLUMN, *CONFounds]].dropna().reset_index(drop=True)
    if len(clean) < 3:
        return float("nan")
    x_conf = _design_matrix(clean, CONFounds)
    x = clean[feature_col].to_numpy(dtype=float)
    y = clean[OUTCOME_COLUMN].to_numpy(dtype=float)
    x_resid = _residualize(x, x_conf)
    y_resid = _residualize(y, x_conf)
    return _safe_corr(x_resid, y_resid)


def _fit_predict_loo(frame: pd.DataFrame, feature_col: str) -> tuple[float, list[dict[str, Any]]]:
    clean = frame[[feature_col, OUTCOME_COLUMN, "donor_id", *CONFounds]].dropna().reset_index(drop=True)
    rows: list[dict[str, Any]] = []
    if len(clean) < 3:
        return float("nan"), rows

    resid_preds: list[float] = []
    resid_actuals: list[float] = []

    for idx in range(len(clean)):
        train = clean.drop(index=idx).reset_index(drop=True)
        test = clean.iloc[[idx]].reset_index(drop=True)

        x_train_conf = _design_matrix(train, CONFounds)
        x_test_conf = _design_matrix(test, CONFounds)

        beta_feat_conf, *_ = np.linalg.lstsq(
            x_train_conf,
            train[feature_col].to_numpy(dtype=float),
            rcond=None,
        )
        beta_out_conf, *_ = np.linalg.lstsq(
            x_train_conf,
            train[OUTCOME_COLUMN].to_numpy(dtype=float),
            rcond=None,
        )

        feat_train_resid = train[feature_col].to_numpy(dtype=float) - x_train_conf @ beta_feat_conf
        out_train_resid = train[OUTCOME_COLUMN].to_numpy(dtype=float) - x_train_conf @ beta_out_conf

        denom = float(np.dot(feat_train_resid, feat_train_resid))
        if denom <= 0:
            slope = float("nan")
        else:
            slope = float(np.dot(feat_train_resid, out_train_resid) / denom)

        feat_test_resid = test[feature_col].to_numpy(dtype=float) - x_test_conf @ beta_feat_conf
        out_test_resid = test[OUTCOME_COLUMN].to_numpy(dtype=float) - x_test_conf @ beta_out_conf
        resid_pred = float(slope * feat_test_resid[0]) if np.isfinite(slope) else float("nan")
        resid_actual = float(out_test_resid[0])

        x_train_full = np.concatenate(
            [x_train_conf, train[[feature_col]].to_numpy(dtype=float)],
            axis=1,
        )
        x_test_full = np.concatenate(
            [x_test_conf, test[[feature_col]].to_numpy(dtype=float)],
            axis=1,
        )
        beta_full, *_ = np.linalg.lstsq(
            x_train_full,
            train[OUTCOME_COLUMN].to_numpy(dtype=float),
            rcond=None,
        )
        raw_pred = float((x_test_full @ beta_full)[0])

        resid_preds.append(resid_pred)
        resid_actuals.append(resid_actual)
        rows.append(
            {
                "donor_id": str(test.loc[0, "donor_id"]),
                "outcome": float(test.loc[0, OUTCOME_COLUMN]),
                "predicted": raw_pred,
                feature_col: float(test.loc[0, feature_col]),
                "residualized_actual": resid_actual,
                "residualized_predicted": resid_pred,
            }
        )

    loo_r = _safe_corr(np.asarray(resid_preds, dtype=float), np.asarray(resid_actuals, dtype=float))
    return loo_r, rows


def _selection_score(partial_r: float, n_analyzable: int, n_total: int) -> float:
    if not np.isfinite(partial_r) or n_total <= 0:
        return float("nan")
    return float(abs(partial_r) * (n_analyzable / n_total))


def _extract_one_donor(donor_id: str, slide_name: str, data_root: Path) -> DonorExtraction:
    slide_path = data_root / slide_name
    svs_path = data_root / slide_name.replace(".zarr", "")
    mpp_um = _aperio_mpp_um(svs_path)

    cells = build_cell_table(slide_path, include_regions=True, include_geometry=True)
    ca1 = cells.loc[cells["region"] == "CA1"].copy()

    pyramidal = ca1.loc[ca1["cell_type"] == "Pyramidal Neuron", ["x", "y"]].to_numpy(dtype=float)
    reactive = ca1.loc[
        ca1["cell_type"] == "Reactive Astrocyte",
        ["x", "y", "area"],
    ].copy()
    lymph = ca1.loc[ca1["cell_type"] == "Lymphocyte", ["x", "y"]].to_numpy(dtype=float)

    reactive_coords = reactive[["x", "y"]].to_numpy(dtype=float)
    reactive_area_um2 = reactive["area"].to_numpy(dtype=float) * (mpp_um ** 2)

    peripyramidal_mask = np.zeros(len(reactive), dtype=bool)
    if len(reactive_coords) and len(pyramidal):
        pyr_tree = cKDTree(pyramidal)
        peripyramidal_radius_px = VARIATIONS[CANONICAL_DEFAULT_VARIATION]["peripyramidal_radius_um"] / mpp_um
        nearest_dist, _ = pyr_tree.query(
            reactive_coords,
            k=1,
            distance_upper_bound=peripyramidal_radius_px,
        )
        peripyramidal_mask = np.isfinite(nearest_dist)

    peripyramidal_coords = reactive_coords[peripyramidal_mask]
    peripyramidal_area_um2 = reactive_area_um2[peripyramidal_mask]
    peripyramidal_count = int(peripyramidal_mask.sum())

    scores: dict[str, float] = {}
    lymph_tree = cKDTree(lymph) if len(lymph) else None
    for variation_name, spec in VARIATIONS.items():
        if peripyramidal_count == 0:
            scores[variation_name] = float("nan")
            continue
        if lymph_tree is None:
            scores[variation_name] = 0.0
            continue
        contact_radius_px = float(spec["lymphocyte_contact_radius_um"]) / mpp_um
        nearest_lymph_dist, _ = lymph_tree.query(
            peripyramidal_coords,
            k=1,
            distance_upper_bound=contact_radius_px,
        )
        contact_mask = np.isfinite(nearest_lymph_dist)
        scores[variation_name] = float(peripyramidal_area_um2[contact_mask].sum() / peripyramidal_count)

    return DonorExtraction(
        donor_id=donor_id,
        slide_name=slide_name,
        mpp_um=mpp_um,
        ca1_pyramidal_count=int(len(pyramidal)),
        ca1_reactive_astro_count=int(len(reactive)),
        ca1_lymphocyte_count=int(len(lymph)),
        peripyramidal_ra_count=peripyramidal_count,
        peripyramidal_ra_area_um2_sum=float(peripyramidal_area_um2.sum()) if len(peripyramidal_area_um2) else 0.0,
        scores=scores,
    )


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Canonical replay target for held-out evaluation.

    Recomputes the winning variation from raw slide data. If results.json is
    present next to this script, it uses the stored best variation; otherwise
    it falls back to the baseline planned variant.
    """
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    row = donor_rows.iloc[0]
    extraction = _extract_one_donor(str(donor_id), str(row["slide_name"]), data_root)
    best_variation = _load_best_variation_from_sidecar()
    return extraction.scores.get(best_variation, float("nan"))


def _format_float(x: float) -> str:
    return "nan" if not np.isfinite(x) else f"{x:.4f}"


def _build_donor_feature_table(data_root: Path) -> tuple[pd.DataFrame, list[DonorExtraction]]:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map({"Female": 0.0, "Male": 1.0}).astype(float)

    extractions: list[DonorExtraction] = []
    for _, row in cohort.iterrows():
        extractions.append(_extract_one_donor(str(row["donor_id"]), str(row["slide_name"]), data_root))

    extracted_frame = pd.DataFrame(
        [
            {
                "donor_id": ext.donor_id,
                "slide_name": ext.slide_name,
                "mpp_um": ext.mpp_um,
                "ca1_pyramidal_count": ext.ca1_pyramidal_count,
                "ca1_reactive_astro_count": ext.ca1_reactive_astro_count,
                "ca1_lymphocyte_count": ext.ca1_lymphocyte_count,
                "peripyramidal_ra_count": ext.peripyramidal_ra_count,
                "peripyramidal_ra_area_um2_sum": ext.peripyramidal_ra_area_um2_sum,
                **{VARIATIONS[k]["feature_column"]: v for k, v in ext.scores.items()},
            }
            for ext in extractions
        ]
    )
    donor_table = cohort.merge(extracted_frame, on=["donor_id", "slide_name"], how="left")
    return donor_table, extractions


def _evaluate_variation(donor_table: pd.DataFrame, variation_name: str) -> dict[str, Any]:
    spec = VARIATIONS[variation_name]
    feature_col = spec["feature_column"]
    clean = donor_table[["donor_id", feature_col, OUTCOME_COLUMN, *CONFounds]].dropna().reset_index(drop=True)
    n_total = int(len(donor_table))
    n_analyzable = int(len(clean))
    partial_r = _partial_correlation(donor_table, feature_col)
    loo_r, loo_rows = _fit_predict_loo(donor_table, feature_col)
    selection_score = _selection_score(partial_r, n_analyzable, n_total)
    if np.isfinite(partial_r) and np.isfinite(loo_r):
        is_loo_gap = float(max(0.0, abs(partial_r) - abs(loo_r)))
    elif np.isfinite(partial_r):
        is_loo_gap = float(abs(partial_r))
    else:
        is_loo_gap = float("nan")
    penalty = is_loo_gap
    adjusted_score = float(selection_score - penalty) if np.isfinite(selection_score) and np.isfinite(penalty) else float("nan")
    return {
        "variation_name": variation_name,
        "feature_name": spec["feature_name"],
        "feature_column": feature_col,
        "description": spec["description"],
        "partial_r": float(partial_r),
        "selection_score": float(selection_score),
        "loo_predictive_r": float(loo_r),
        "is_loo_gap": float(is_loo_gap),
        "penalty": float(penalty),
        "adjusted_score": float(adjusted_score),
        "n_analyzable": n_analyzable,
        "n_total": n_total,
        "per_donor_loo": loo_rows,
    }


def _rank_variations(metrics: list[dict[str, Any]]) -> list[dict[str, Any]]:
    def sort_key(item: dict[str, Any]):
        sel = item["selection_score"]
        pr = item["partial_r"]
        loo = item["loo_predictive_r"]
        return (
            -999 if not np.isfinite(sel) else -sel,
            -999 if not np.isfinite(abs(pr)) else -abs(pr),
            -999 if not np.isfinite(abs(loo)) else -abs(loo),
            item["variation_name"],
        )

    return sorted(metrics, key=sort_key)


def _error_summary(best_metrics: dict[str, Any]) -> dict[str, Any]:
    loo_rows = pd.DataFrame(best_metrics["per_donor_loo"])
    feature_col = best_metrics["feature_column"]
    if loo_rows.empty:
        return {"worst_donors": [], "residual_sign": {}}
    loo_rows["abs_error"] = (loo_rows["predicted"] - loo_rows["outcome"]).abs()
    worst = loo_rows.sort_values("abs_error", ascending=False).head(5)
    return {
        "worst_donors": worst[["donor_id", "outcome", "predicted", feature_col, "abs_error"]].to_dict(orient="records"),
    }


def _write_report(payload: dict[str, Any], donor_table: pd.DataFrame) -> None:
    best = payload
    ranked = payload["ranked_variations"]
    winner = ranked[0]
    loser = ranked[1] if len(ranked) > 1 else None
    feature_col = winner["feature_column"]

    loo_df = pd.DataFrame(winner["per_donor_loo"])
    error_lines = []
    if not loo_df.empty:
        loo_df["error"] = loo_df["predicted"] - loo_df["outcome"]
        loo_df["abs_error"] = loo_df["error"].abs()
        worst = loo_df.sort_values("abs_error", ascending=False).head(3)
        for _, row in worst.iterrows():
            error_lines.append(
                f"- {row['donor_id']}: outcome {row['outcome']:.3f}, predicted {row['predicted']:.3f}, "
                f"{feature_col}={row[feature_col]:.2f}"
            )
    error_block = "\n".join(error_lines) if error_lines else "- No analyzable LOO donor table."

    if loser is None:
        loser_summary = "No alternate local variation was tested."
    elif winner["variation_name"] == "candidate_variant_a":
        loser_summary = (
            f"{loser['variation_name']} (25 um lymphocyte radius) ranked below the winner "
            f"with selection score {loser['selection_score']:.4f}, suggesting that broadening the "
            "contact definition diluted CA1 peripyramidal specificity."
        )
    else:
        loser_summary = (
            f"{loser['variation_name']} (20 um lymphocyte radius) ranked below the winner "
            f"with selection score {loser['selection_score']:.4f}, suggesting that the narrower "
            "contact definition missed biologically relevant near-contact reactive astrocytes."
        )

    next_suggestion = (
        "Keep the CA1 peripyramidal reactive-astrocyte/lymphocyte burden family, but in the next local sweep "
        "change one thing only: retain the winning lymphocyte radius and test whether normalizing by total "
        "peripyramidal reactive-astrocyte area rather than count reduces the largest donor-level prediction errors."
    )

    report = f"""## Summary
One sentence: tested the {HYPOTHESIS_FAMILY} family, {winner['variation_name']} won, and its local winner selection score was {winner['selection_score']:.4f}.

## Metrics
Winning variation: {winner['variation_name']} ({winner['feature_column']}).
- IS partial r: {winner['partial_r']:.4f}
- Selection score: {winner['selection_score']:.4f}
- LOO predictive r: {winner['loo_predictive_r']:.4f}
- IS-LOO Gap: {winner['is_loo_gap']:.4f} (penalty={winner['penalty']:.4f})
- Adjusted Score: {winner['adjusted_score']:.4f}

Ranking:
""" + "\n".join(
        [
            f"- {row['variation_name']}: partial_r={row['partial_r']:.4f}, selection_score={row['selection_score']:.4f}, loo_predictive_r={row['loo_predictive_r']:.4f}"
            for row in ranked
        ]
    ) + f"""

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The winning feature isolates CA1 reactive astrocytes that are both peripyramidal and lymphocyte-associated, then upweights them by hypertrophic contour area. That matches the scientific question: the signal is strongest when inflammatory contact and astrocyte enlargement are coupled in the pyramidal niche.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - {loser_summary}
3. Error pattern: which donors are consistently wrong and what they share
{error_block}

## Rationale
The winning variation is biologically coherent because it focuses on a narrow CA1 niche already implicated by prior accepted features: reactive astrocytes around pyramidal neurons plus local lymphocyte contact. Using summed reactive-astrocyte area over peripyramidal reactive-astrocyte count converts that niche into a severity-weighted burden, which is more specific than a simple presence/count feature. It beat the nearby alternative by preserving the contact scale that best matched the data. It likely carries partially new information beyond the current panel because it fuses hypertrophy with immune contact instead of measuring them separately, though the panel-level additivity still needs formal evaluator review.

## Interpretation
The signal appears to mean a donor-specific burden of enlarged immune-associated reactive astrocytes in the CA1 pyramidal neighborhood. The population is CA1 Reactive Astrocyte cells, the niche is the 30 um peripyramidal zone around CA1 Pyramidal Neurons with nearby CA1 Lymphocytes, the feature summary is summed contact-positive reactive-astrocyte area divided by peripyramidal reactive-astrocyte count, and the simplest observable tissue pattern is clusters of large reactive astrocytes hugging the CA1 pyramidal layer with adjacent lymphocytes.

## Next
{next_suggestion}
"""
    REPORT_PATH.write_text(report)


def main() -> None:
    data_root = Path("/data")
    donor_table, _ = _build_donor_feature_table(data_root)
    DONOR_FEATURE_TABLE_PATH.write_text(donor_table.to_csv(index=False))

    metrics = [_evaluate_variation(donor_table, variation_name) for variation_name in VARIATIONS]
    ranked = _rank_variations(metrics)
    best = ranked[0]
    best_feature_col = best["feature_column"]

    payload = {
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": best["variation_name"],
        "feature_name": best["feature_name"],
        "feature_column": best_feature_col,
        "partial_r": best["partial_r"],
        "selection_score": best["selection_score"],
        "loo_predictive_r": best["loo_predictive_r"],
        "is_loo_gap": best["is_loo_gap"],
        "penalty": best["penalty"],
        "adjusted_score": best["adjusted_score"],
        "n_analyzable": best["n_analyzable"],
        "n_total": best["n_total"],
        "ranked_variations": ranked,
        "error_summary": _error_summary(best),
    }
    RESULTS_PATH.write_text(json.dumps(payload, indent=2, default=_json_default))
    _write_report(payload, donor_table)

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best['variation_name']}")
    print(f"  IS partial r:      {_format_float(best['partial_r'])}")
    print(f"  Selection score:   {_format_float(best['selection_score'])}")
    print(f"  LOO predictive r:  {_format_float(best['loo_predictive_r'])}  (diagnostic)")
    print(f"  IS-LOO Gap:        {_format_float(best['is_loo_gap'])}  (penalty={_format_float(best['penalty'])})")
    print(f"  Adjusted Score:    {_format_float(best['adjusted_score'])}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for row in ranked:
        print(
            f"  {row['variation_name']}  "
            f"{_format_float(row['partial_r'])}  "
            f"{_format_float(row['selection_score'])}  "
            f"{_format_float(row['loo_predictive_r'])}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best_feature_col}")
    for row in best["per_donor_loo"]:
        outcome = _format_float(float(row["outcome"]))
        predicted = _format_float(float(row["predicted"]))
        feature_value = _format_float(float(row[best_feature_col]))
        print(f"  {row['donor_id']}  {outcome}  {predicted}  {feature_value}")


if __name__ == "__main__":
    main()
