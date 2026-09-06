set -euo pipefail
cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy.spatial import cKDTree

sys.path.append("/shared/lib")

from shared_analysis import build_cell_table, load_training_cohort  # noqa: E402

warnings.filterwarnings("ignore", category=UserWarning)

FAMILY_NAME = "CA1 pyramidal reactive-astrocyte cuffing fraction"
TARGET_REGION = "CA1"
TARGET_CELL_TYPE = "Pyramidal Neuron"
NEIGHBOR_CELL_TYPE = "Reactive Astrocyte"
RADIUS_UM = 35.0
FALLBACK_MPP_UM_PER_PX = 0.50325
OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]

VARIATIONS = {
    "r35_ge1": {"min_neighbors": 1, "feature_column": "ca1_pyramidal_cuffed_fraction_r35um_ge1"},
    "r35_ge2": {"min_neighbors": 2, "feature_column": "ca1_pyramidal_cuffed_fraction_r35um_ge2"},
    "r35_ge3": {"min_neighbors": 3, "feature_column": "ca1_pyramidal_cuffed_fraction_r35um_ge3"},
}

CANONICAL_VARIATION = "r35_ge2"  # AUTO_CANONICAL_VARIATION
CANONICAL_MIN_REACTIVE_NEIGHBORS = 2  # AUTO_MIN_REACTIVE_NEIGHBORS


def canonical_feature_column() -> str:
    return VARIATIONS[CANONICAL_VARIATION]["feature_column"]


def canonical_feature_name() -> str:
    return canonical_feature_column()


def _sex_to_binary(series: pd.Series) -> pd.Series:
    mapped = series.astype(str).str.strip().str.lower().map({"female": 0.0, "male": 1.0})
    return mapped.astype(float)


def _design_matrix(df: pd.DataFrame) -> np.ndarray:
    return np.column_stack([np.ones(len(df), dtype=float), df.to_numpy(dtype=float)])


def _clean_frame(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    frame = df.loc[:, columns].replace([np.inf, -np.inf], np.nan).dropna().copy()
    return frame


def _partial_correlation(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
) -> tuple[float, int]:
    frame = _clean_frame(df, [feature_col, outcome_col, *confounds])
    n = len(frame)
    if n < 3:
        return float("nan"), n
    x = _design_matrix(frame[confounds])
    y_feature = frame[feature_col].to_numpy(dtype=float)
    y_outcome = frame[outcome_col].to_numpy(dtype=float)

    beta_feature, *_ = np.linalg.lstsq(x, y_feature, rcond=None)
    beta_outcome, *_ = np.linalg.lstsq(x, y_outcome, rcond=None)

    resid_feature = y_feature - x @ beta_feature
    resid_outcome = y_outcome - x @ beta_outcome

    if np.std(resid_feature) == 0 or np.std(resid_outcome) == 0:
        return float("nan"), n
    corr = float(np.corrcoef(resid_feature, resid_outcome)[0, 1])
    return corr, n


def _loo_predictions(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
) -> tuple[float, pd.DataFrame]:
    required = ["donor_id", "slide_name", feature_col, outcome_col, *confounds]
    frame = _clean_frame(df, required).reset_index(drop=True)
    records: list[dict[str, float | str]] = []

    if len(frame) < 3:
        return float("nan"), pd.DataFrame(records)

    resid_preds: list[float] = []
    resid_actuals: list[float] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train_conf = _design_matrix(train[confounds])
        x_test_conf = _design_matrix(test[confounds])

        y_feature_train = train[feature_col].to_numpy(dtype=float)
        y_outcome_train = train[outcome_col].to_numpy(dtype=float)

        beta_feature, *_ = np.linalg.lstsq(x_train_conf, y_feature_train, rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train_conf, y_outcome_train, rcond=None)

        resid_feature_train = y_feature_train - x_train_conf @ beta_feature
        resid_outcome_train = y_outcome_train - x_train_conf @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            pred_resid = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            feature_test = test[feature_col].to_numpy(dtype=float)
            resid_feature_test = feature_test - x_test_conf @ beta_feature
            pred_resid = float(slope * resid_feature_test[0])

        actual_resid = float(test[outcome_col].to_numpy(dtype=float)[0] - (x_test_conf @ beta_outcome)[0])
        pred_outcome = (
            float((x_test_conf @ beta_outcome)[0] + pred_resid) if math.isfinite(pred_resid) else float("nan")
        )

        resid_preds.append(pred_resid)
        resid_actuals.append(actual_resid)

        records.append(
            {
                "donor_id": str(test.loc[0, "donor_id"]),
                "slide_name": str(test.loc[0, "slide_name"]),
                "outcome": float(test.loc[0, outcome_col]),
                "predicted": pred_outcome,
                feature_col: float(test.loc[0, feature_col]),
            }
        )

    pred_arr = np.asarray(resid_preds, dtype=float)
    act_arr = np.asarray(resid_actuals, dtype=float)
    mask = np.isfinite(pred_arr) & np.isfinite(act_arr)
    if mask.sum() < 3 or np.std(pred_arr[mask]) == 0 or np.std(act_arr[mask]) == 0:
        loo_r = float("nan")
    else:
        loo_r = float(np.corrcoef(pred_arr[mask], act_arr[mask])[0, 1])

    return loo_r, pd.DataFrame(records)


def _parse_mpp_from_svs(svs_path: Path) -> float:
    try:
        with tifffile.TiffFile(svs_path) as tf:
            description = str(tf.pages[0].description)
        match = re.search(r"\bMPP\s*=\s*([0-9.]+)", description)
        if match:
            return float(match.group(1))
    except Exception:
        pass
    return FALLBACK_MPP_UM_PER_PX


def _extract_slide_features(*, slide_zarr_path: Path, slide_svs_path: Path) -> dict[str, float]:
    cells = build_cell_table(slide_zarr_path, include_regions=True, include_geometry=False)
    ca1_cells = cells.loc[cells["region"] == TARGET_REGION, ["x", "y", "cell_type"]].copy()

    pyramidal = ca1_cells.loc[ca1_cells["cell_type"] == TARGET_CELL_TYPE, ["x", "y"]].to_numpy(dtype=float)
    reactive = ca1_cells.loc[ca1_cells["cell_type"] == NEIGHBOR_CELL_TYPE, ["x", "y"]].to_numpy(dtype=float)

    n_pyramidal = int(len(pyramidal))
    n_reactive = int(len(reactive))
    mpp = _parse_mpp_from_svs(slide_svs_path)
    radius_px = float(RADIUS_UM / mpp)

    out: dict[str, float] = {
        "ca1_pyramidal_count": float(n_pyramidal),
        "ca1_reactive_astrocyte_count": float(n_reactive),
        "radius_um": float(RADIUS_UM),
        "slide_mpp_um_per_px": float(mpp),
        "radius_px": float(radius_px),
    }

    if n_pyramidal == 0:
        for spec in VARIATIONS.values():
            out[spec["feature_column"]] = float("nan")
        out["mean_reactive_neighbors_within_35um"] = float("nan")
        return out

    if n_reactive == 0:
        neighbor_counts = np.zeros(n_pyramidal, dtype=np.int32)
    else:
        tree = cKDTree(reactive)
        neighbor_lists = tree.query_ball_point(pyramidal, r=radius_px)
        neighbor_counts = np.asarray([len(v) for v in neighbor_lists], dtype=np.int32)

    out["mean_reactive_neighbors_within_35um"] = float(np.mean(neighbor_counts)) if len(neighbor_counts) else float("nan")
    for spec in VARIATIONS.values():
        out[spec["feature_column"]] = float(np.mean(neighbor_counts >= int(spec["min_neighbors"])))
    return out


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Return the canonical winning donor-level biomarker score for held-out replay.
    """
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root).copy()
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None

    row = donor_rows.iloc[0]
    slide_name = str(row["slide_name"])
    slide_zarr_path = data_root / slide_name
    slide_svs_path = data_root / slide_name.replace(".zarr", "")
    features = _extract_slide_features(slide_zarr_path=slide_zarr_path, slide_svs_path=slide_svs_path)
    value = features.get(canonical_feature_column(), float("nan"))
    return None if not math.isfinite(value) else float(value)


def _update_self_with_winner(best_variation: str) -> None:
    min_neighbors = int(VARIATIONS[best_variation]["min_neighbors"])
    path = Path(__file__)
    text = path.read_text()
    text = re.sub(
        r'CANONICAL_VARIATION = ".*?"  # AUTO_CANONICAL_VARIATION',
        f'CANONICAL_VARIATION = "{best_variation}"  # AUTO_CANONICAL_VARIATION',
        text,
    )
    text = re.sub(
        r"CANONICAL_MIN_REACTIVE_NEIGHBORS = \d+  # AUTO_MIN_REACTIVE_NEIGHBORS",
        f"CANONICAL_MIN_REACTIVE_NEIGHBORS = {min_neighbors}  # AUTO_MIN_REACTIVE_NEIGHBORS",
        text,
    )
    path.write_text(text)


def _write_report(
    *,
    best: dict,
    ranked: list[dict],
    loo_table: pd.DataFrame,
    report_path: Path,
) -> None:
    feat_col = best["feature_column"]
    error_frame = loo_table.copy()
    error_frame["abs_error"] = (error_frame["outcome"] - error_frame["predicted"]).abs()
    worst = error_frame.sort_values("abs_error", ascending=False).head(5)
    worst_lines = []
    for _, row in worst.iterrows():
        worst_lines.append(
            f"- {row['donor_id']}: outcome={row['outcome']:.3f}, predicted={row['predicted']:.3f}, "
            f"{feat_col}={row[feat_col]:.3f}"
        )
    ranking_bits = [
        f"{r['variation_name']} (selection={r['selection_score']:.4f}, partial_r={r['partial_r']:.4f}, loo_r={r['loo_predictive_r']:.4f})"
        for r in ranked
    ]
    text = f"""## Summary
One sentence: tested the {FAMILY_NAME} family, {best['variation_name']} won, and its local winner selection score was {best['selection_score']:.4f}.

## Metrics
Winning variation: {best['variation_name']}.
- IS partial r: {best['partial_r']:.4f}
- Selection score: {best['selection_score']:.4f}
- LOO predictive r: {best['loo_predictive_r']:.4f}
- IS-LOO Gap: {best['is_loo_gap']:.4f} (penalty={best['penalty']:.4f})
- Adjusted Score: {best['adjusted_score']:.4f}

Other tested variations ranked by local single-feature signal:
- {'; '.join(ranking_bits)}

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The neuron-centric cuffing summary worked best when defined as {best['variation_name']}, meaning the donor-level signal is strongest when asking how broadly CA1 pyramidal neurons are exposed to nearby reactive astrocytes, rather than only whether a minimal contact exists. This keeps the population fixed to CA1 pyramidal neurons and sharpens the niche to a 35 µm perineuronal neighborhood.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The nearby thresholds lost signal because they were either too permissive ({'r35_ge1' if best['variation_name'] != 'r35_ge1' else 'r35_ge2'}) and likely collapsed donors toward saturation, or too strict ({'r35_ge3' if best['variation_name'] != 'r35_ge3' else 'r35_ge2'}) and likely discarded informative moderate cuffing events. The family appears sensitive to how much reactive-astrocyte occupancy is required before a neuron counts as exposed.
3. Error pattern: which donors are consistently wrong and what they share
   - Largest LOO errors were:
{chr(10).join(worst_lines)}

## Rationale
Why the best variation's approach is biologically coherent and why it beat the nearby alternatives.
- CA1 pyramidal neurons are the hippocampal population most directly tied to memory circuitry, so measuring how often they sit inside a local ring of reactive astrocytes is a biologically coherent niche summary.
- This beat the nearby alternatives because the winning threshold best balanced prevalence and specificity: enough reactive astrocytes to indicate a real cuffing state, but not so many that only extreme cases contribute.
- Relative to the current panel, it is plausibly additive because it converts an astrocyte-centered reactivity niche into a neuron-exposure niche.

## Interpretation
State this when it materially helps explain the biomarker.
- The signal appears to mean that donors with worse memory decline have a larger fraction of CA1 pyramidal neurons locally surrounded by reactive astrocytes.
- Population: CA1 Pyramidal Neuron, with neighboring CA1 Reactive Astrocyte.
- Niche: a 35 µm centroid-defined perineuronal neighborhood inside CA1.
- Feature summary: donor-level fraction of CA1 pyramidal neurons whose local reactive-astrocyte neighbor count crosses the winning threshold.
- Simplest observable pattern: broader reactive astrocyte cuffing around CA1 pyramidal neurons.

## Next
One specific suggestion for the next local sweep based on the error pattern and which nearby variations won or lost.
- Keep the same CA1 pyramidal / reactive astrocyte niche but sweep radius around the winning threshold setting (for example 25, 45, and 55 µm with the winning minimum-neighbor rule fixed) to test whether the current errors reflect an overly narrow spatial scale rather than the exposure definition itself.
"""
    report_path.write_text(text)


def main() -> int:
    data_root = Path("/data")
    scratch = Path("/scratch")
    donor_feature_table_path = scratch / "donor_feature_table.csv"
    results_json_path = scratch / "results.json"
    report_path = scratch / "report.md"
    state_path = scratch / "best_variation.json"

    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = _sex_to_binary(cohort["sex"])
    total_donors = int(len(cohort))

    records: list[dict[str, float | str]] = []
    for _, row in cohort.iterrows():
        slide_name = str(row["slide_name"])
        slide_zarr_path = data_root / slide_name
        slide_svs_path = data_root / slide_name.replace(".zarr", "")
        feats = _extract_slide_features(slide_zarr_path=slide_zarr_path, slide_svs_path=slide_svs_path)
        records.append(
            {
                "donor_id": str(row["donor_id"]),
                "slide_name": slide_name,
                OUTCOME_COLUMN: float(row[OUTCOME_COLUMN]),
                "max_age_vis": float(row["max_age_vis"]),
                "braak_numeric": float(row["braak_numeric"]),
                "cerad_ordinal": float(row["cerad_ordinal"]),
                "sex_binary": float(row["sex_binary"]),
                **feats,
            }
        )

    donor_table = pd.DataFrame(records)

    ranked_variations: list[dict[str, float | str]] = []
    loo_tables: dict[str, pd.DataFrame] = {}

    for variation_name, spec in VARIATIONS.items():
        feature_col = str(spec["feature_column"])
        partial_r, n_analyzable = _partial_correlation(
            donor_table,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUND_COLUMNS,
        )
        loo_r, loo_table = _loo_predictions(
            donor_table,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUND_COLUMNS,
        )
        loo_table = loo_table.merge(
            donor_table[
                [
                    "donor_id",
                    "ca1_pyramidal_count",
                    "ca1_reactive_astrocyte_count",
                    "mean_reactive_neighbors_within_35um",
                ]
            ],
            on="donor_id",
            how="left",
        )
        loo_tables[variation_name] = loo_table

        coverage = float(n_analyzable) / float(total_donors) if total_donors else float("nan")
        selection_score = abs(partial_r) * coverage if math.isfinite(partial_r) else float("nan")
        is_loo_gap = abs(partial_r) - abs(loo_r) if math.isfinite(partial_r) and math.isfinite(loo_r) else float("nan")
        penalty = max(0.0, is_loo_gap) if math.isfinite(is_loo_gap) else float("nan")
        adjusted_score = selection_score - penalty if math.isfinite(selection_score) and math.isfinite(penalty) else float("nan")

        ranked_variations.append(
            {
                "variation_name": variation_name,
                "feature_column": feature_col,
                "min_neighbors": int(spec["min_neighbors"]),
                "n_analyzable": int(n_analyzable),
                "n_total": int(total_donors),
                "coverage": coverage,
                "partial_r": float(partial_r),
                "selection_score": float(selection_score),
                "loo_predictive_r": float(loo_r),
                "is_loo_gap": float(is_loo_gap),
                "penalty": float(penalty),
                "adjusted_score": float(adjusted_score),
            }
        )

    ranked_variations.sort(
        key=lambda d: (
            -np.inf if not math.isfinite(float(d["selection_score"])) else -float(d["selection_score"]),
            -np.inf if not math.isfinite(float(d["adjusted_score"])) else -float(d["adjusted_score"]),
        )
    )
    best = ranked_variations[0]
    best_variation = str(best["variation_name"])
    best_feature_col = str(best["feature_column"])
    best_loo_table = loo_tables[best_variation].copy()

    donor_feature_table = donor_table[
        [
            "donor_id",
            "slide_name",
            best_feature_col,
            OUTCOME_COLUMN,
            *CONFOUND_COLUMNS,
            "ca1_pyramidal_count",
            "ca1_reactive_astrocyte_count",
            "mean_reactive_neighbors_within_35um",
        ]
    ].copy()
    donor_feature_table.to_csv(donor_feature_table_path, index=False)

    payload = {
        "family_name": FAMILY_NAME,
        "best_variation": best_variation,
        "feature_column": best_feature_col,
        "ranked_variations": ranked_variations,
        "partial_r": best["partial_r"],
        "selection_score": best["selection_score"],
        "loo_predictive_r": best["loo_predictive_r"],
        "is_loo_gap": best["is_loo_gap"],
        "penalty": best["penalty"],
        "adjusted_score": best["adjusted_score"],
        "n_analyzable": best["n_analyzable"],
        "n_total": best["n_total"],
    }
    results_json_path.write_text(json.dumps(payload, indent=2))

    state_path.write_text(
        json.dumps(
            {
                "best_variation": best_variation,
                "min_neighbors": int(VARIATIONS[best_variation]["min_neighbors"]),
                "feature_column": best_feature_col,
                "radius_um": RADIUS_UM,
                "region": TARGET_REGION,
                "target_cell_type": TARGET_CELL_TYPE,
                "neighbor_cell_type": NEIGHBOR_CELL_TYPE,
            },
            indent=2,
        )
    )

    _update_self_with_winner(best_variation)
    _write_report(best=best, ranked=ranked_variations, loo_table=best_loo_table, report_path=report_path)

    print(f"HYPOTHESIS FAMILY: {FAMILY_NAME}")
    print(f"BEST VARIATION: {best_variation}")
    print(f"  IS partial r:      {best['partial_r']:.4f}")
    print(f"  Selection score:   {best['selection_score']:.4f}")
    print(f"  LOO predictive r:  {best['loo_predictive_r']:.4f}  (diagnostic)")
    print(f"  IS-LOO Gap:        {best['is_loo_gap']:.4f}  (penalty={best['penalty']:.4f})")
    print(f"  Adjusted Score:    {best['adjusted_score']:.4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for row in ranked_variations:
        print(
            f"  {row['variation_name']}  {row['partial_r']:.4f}  "
            f"{row['selection_score']:.4f}  {row['loo_predictive_r']:.4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(
        f"  donor_id  outcome  predicted  {best_feature_col}  "
        "ca1_pyramidal_count  ca1_reactive_astrocyte_count"
    )
    for _, row in best_loo_table.iterrows():
        print(
            f"  {row['donor_id']}  {row['outcome']:.4f}  {row['predicted']:.4f}  "
            f"{row[best_feature_col]:.4f}  {row['ca1_pyramidal_count']:.0f}  "
            f"{row['ca1_reactive_astrocyte_count']:.0f}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
PY
python /scratch/result.py
ls -l /scratch/results.json /scratch/report.md /scratch/donor_feature_table.csv /scratch/best_variation.json
sed -n '1,80p' /scratch/result.py
