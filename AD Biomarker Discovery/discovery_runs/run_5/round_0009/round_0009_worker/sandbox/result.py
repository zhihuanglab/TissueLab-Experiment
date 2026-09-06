from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore", message=".*not recognized as a component of a Zarr hierarchy.*")

# Make shared helpers importable in the sandbox runtime.
sys.path.insert(0, "/shared/lib")

from shared_analysis import build_cell_table, load_training_cohort  # type: ignore
from shared_analysis.stats import (  # type: ignore
    partial_correlation,
    residualized_loo_predictive_correlation,
)

HYPOTHESIS_FAMILY = "CA1 reactive-astro-exposed pyramidal lymphocyte enrichment within 25um"

# Default canonical replay target; if results.json is present next to result.py,
# compute_donor_score will resolve the actual winning variation from that file.
CANONICAL_VARIATION = "candidate_variant_a"

VARIATIONS: dict[str, dict[str, Any]] = {
    "candidate_variant_a": {
        "description": "Fraction of CA1 reactive-astro-exposed pyramidal neurons with at least 1 lymphocyte within 25 um.",
        "lymphocyte_threshold": 1,
        "feature_column": "ca1_ra_exposed_pyramidal_lym_ge1_fraction_25um",
    },
    "candidate_variant_b": {
        "description": "Fraction of CA1 reactive-astro-exposed pyramidal neurons with at least 2 lymphocytes within 25 um.",
        "lymphocyte_threshold": 2,
        "feature_column": "ca1_ra_exposed_pyramidal_lym_ge2_fraction_25um",
    },
}

FEATURE_COLUMN = VARIATIONS[CANONICAL_VARIATION]["feature_column"]
OUTCOME_COLUMN = "slope_zmem0"
CONFOUNDS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
RADIUS_UM = 25.0
MPP_FALLBACK = 0.5032


def _safe_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        val = float(value)
    except Exception:
        return None
    return val if math.isfinite(val) else None


def _json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _json_ready(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_json_ready(v) for v in value]
    if isinstance(value, tuple):
        return [_json_ready(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        val = float(value)
        return val if math.isfinite(val) else None
    if pd.isna(value):
        return None
    return value


def _resolve_results_path() -> Path:
    here = Path(__file__).resolve().parent
    scratch_path = here / "results.json"
    return scratch_path


def resolve_canonical_variation() -> str:
    results_path = _resolve_results_path()
    if results_path.exists():
        try:
            payload = json.loads(results_path.read_text())
            best = payload.get("best_variation")
            if best in VARIATIONS:
                return str(best)
        except Exception:
            pass
    return CANONICAL_VARIATION


def read_slide_mpp(svs_path: Path) -> float:
    try:
        import openslide  # type: ignore

        slide = openslide.OpenSlide(str(svs_path))
        raw = slide.properties.get("aperio.MPP") or slide.properties.get(openslide.PROPERTY_NAME_MPP_X)
        if raw is not None:
            mpp = float(raw)
            if math.isfinite(mpp) and mpp > 0:
                return mpp
    except Exception:
        pass
    return MPP_FALLBACK


def _count_neighbors(tree: cKDTree | None, query_points: np.ndarray, radius_px: float) -> np.ndarray:
    query_points = np.asarray(query_points, dtype=np.float32)
    if tree is None or len(query_points) == 0:
        return np.zeros(len(query_points), dtype=np.int32)
    try:
        counts = tree.query_ball_point(query_points, r=float(radius_px), return_length=True)
        return np.asarray(counts, dtype=np.int32)
    except TypeError:
        hits = tree.query_ball_point(query_points, r=float(radius_px))
        return np.asarray([len(x) for x in hits], dtype=np.int32)


def compute_slide_family_features(*, slide_path: Path) -> dict[str, Any]:
    svs_path = slide_path.with_suffix("")
    mpp = read_slide_mpp(svs_path)
    radius_px = float(RADIUS_UM / mpp)

    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    ca1 = cells.loc[cells["region"] == "CA1", ["x", "y", "cell_type"]].copy()

    row: dict[str, Any] = {
        "slide_name": slide_path.name,
        "mpp_um_per_px": float(mpp),
        "radius_px_25um": float(radius_px),
        "n_cells_ca1": int(len(ca1)),
    }

    if ca1.empty:
        row.update(
            {
                "n_pyramidal_ca1": 0,
                "n_reactive_astro_ca1": 0,
                "n_lymphocyte_ca1": 0,
                "n_ra_exposed_pyramidal_ca1": 0,
                "median_lymphocytes_per_exposed_pyramidal_25um": np.nan,
                "ca1_missing": 1,
            }
        )
        for spec in VARIATIONS.values():
            row[spec["feature_column"]] = np.nan
        return row

    pyramidal = ca1.loc[ca1["cell_type"] == "Pyramidal Neuron", ["x", "y"]].to_numpy(dtype=np.float32)
    reactive_astro = ca1.loc[ca1["cell_type"] == "Reactive Astrocyte", ["x", "y"]].to_numpy(dtype=np.float32)
    lymphocytes = ca1.loc[ca1["cell_type"] == "Lymphocyte", ["x", "y"]].to_numpy(dtype=np.float32)

    row["n_pyramidal_ca1"] = int(len(pyramidal))
    row["n_reactive_astro_ca1"] = int(len(reactive_astro))
    row["n_lymphocyte_ca1"] = int(len(lymphocytes))
    row["ca1_missing"] = 0

    if len(pyramidal) == 0 or len(reactive_astro) == 0:
        row["n_ra_exposed_pyramidal_ca1"] = 0
        row["median_lymphocytes_per_exposed_pyramidal_25um"] = np.nan
        for spec in VARIATIONS.values():
            row[spec["feature_column"]] = np.nan
        return row

    reactive_tree = cKDTree(reactive_astro)
    exposure_counts = _count_neighbors(reactive_tree, pyramidal, radius_px)
    exposed_mask = exposure_counts > 0
    n_exposed = int(exposed_mask.sum())
    row["n_ra_exposed_pyramidal_ca1"] = n_exposed

    if n_exposed == 0:
        row["median_lymphocytes_per_exposed_pyramidal_25um"] = np.nan
        for spec in VARIATIONS.values():
            row[spec["feature_column"]] = np.nan
        return row

    exposed_pyramidal = pyramidal[exposed_mask]
    lymph_tree = cKDTree(lymphocytes) if len(lymphocytes) else None
    lymph_counts = _count_neighbors(lymph_tree, exposed_pyramidal, radius_px)
    row["median_lymphocytes_per_exposed_pyramidal_25um"] = float(np.median(lymph_counts)) if len(lymph_counts) else np.nan

    for spec in VARIATIONS.values():
        threshold = int(spec["lymphocyte_threshold"])
        row[spec["feature_column"]] = float(np.mean(lymph_counts >= threshold))

    return row


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Return the winning biomarker value for one donor.

    Replays from raw zarr + paired SVS metadata under data_root and resolves the
    canonical winning variation from results.json next to this script when
    available.
    """
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    feature_row = compute_slide_family_features(slide_path=data_root / slide_name)
    variation = resolve_canonical_variation()
    feature_column = VARIATIONS[variation]["feature_column"]
    score = feature_row.get(feature_column)
    return _safe_float(score)


def build_loo_predictions(df: pd.DataFrame, *, feature_col: str) -> pd.DataFrame:
    frame = df.loc[:, ["donor_id", OUTCOME_COLUMN, feature_col, "n_ra_exposed_pyramidal_ca1", *CONFOUNDS]].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    rows: list[dict[str, Any]] = []
    if frame.empty:
        return pd.DataFrame(columns=["donor_id", "outcome", "predicted", feature_col, "n_ra_exposed_pyramidal_ca1"])

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train = np.column_stack([np.ones(len(train)), train[CONFOUNDS].astype(float).to_numpy()])
        x_test = np.column_stack([np.ones(len(test)), test[CONFOUNDS].astype(float).to_numpy()])

        y_feature_train = train[feature_col].to_numpy(dtype=float)
        y_outcome_train = train[OUTCOME_COLUMN].to_numpy(dtype=float)

        beta_feature, *_ = np.linalg.lstsq(x_train, y_feature_train, rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train, y_outcome_train, rcond=None)

        resid_feature_train = y_feature_train - x_train @ beta_feature
        resid_outcome_train = y_outcome_train - x_train @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            pred_outcome = np.nan
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            resid_feature_test = test[feature_col].to_numpy(dtype=float) - x_test @ beta_feature
            pred_resid = float(slope * resid_feature_test[0])
            pred_outcome = float((x_test @ beta_outcome)[0] + pred_resid)

        rows.append(
            {
                "donor_id": str(test.loc[0, "donor_id"]),
                "outcome": float(test.loc[0, OUTCOME_COLUMN]),
                "predicted": pred_outcome,
                feature_col: float(test.loc[0, feature_col]),
                "n_ra_exposed_pyramidal_ca1": int(test.loc[0, "n_ra_exposed_pyramidal_ca1"]),
            }
        )
    return pd.DataFrame(rows)


def evaluate_variation(df: pd.DataFrame, *, variation_name: str) -> dict[str, Any]:
    spec = VARIATIONS[variation_name]
    feature_col = spec["feature_column"]

    usable = df.loc[:, ["donor_id", OUTCOME_COLUMN, feature_col, *CONFOUNDS]].replace([np.inf, -np.inf], np.nan).dropna()
    n_total = int(len(df))
    n_analyzable = int(len(usable))

    metrics = partial_correlation(
        df,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
    )
    partial_r = float(metrics["partial_r"]) if metrics["partial_r"] == metrics["partial_r"] else np.nan
    loo_predictive_r = float(
        residualized_loo_predictive_correlation(
            df,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUNDS,
        )
    ) if n_analyzable >= 6 else np.nan

    coverage = float(n_analyzable / n_total) if n_total else np.nan
    selection_score = abs(partial_r) * coverage if math.isfinite(partial_r) and math.isfinite(coverage) else np.nan
    is_loo_gap = (abs(partial_r) - abs(loo_predictive_r)) if (math.isfinite(partial_r) and math.isfinite(loo_predictive_r)) else np.nan
    penalty = max(0.0, float(is_loo_gap)) if math.isfinite(is_loo_gap) else 0.0
    adjusted_score = (selection_score - penalty) if math.isfinite(selection_score) else np.nan

    loo_table = build_loo_predictions(df, feature_col=feature_col)

    return {
        "variation_name": variation_name,
        "description": spec["description"],
        "feature_column": feature_col,
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "coverage_fraction": coverage,
        "partial_r": partial_r,
        "p_value": float(metrics["p_value"]) if metrics["p_value"] == metrics["p_value"] else np.nan,
        "selection_score": selection_score,
        "loo_predictive_r": loo_predictive_r,
        "is_loo_gap": is_loo_gap,
        "penalty": penalty,
        "adjusted_score": adjusted_score,
        "loo_table": loo_table,
    }


def main() -> None:
    data_root = Path("/data")
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map({"Female": 0.0, "Male": 1.0})

    feature_rows: list[dict[str, Any]] = []
    for _, row in cohort.iterrows():
        donor_id = str(row["donor_id"])
        slide_name = str(row["slide_name"])
        feature_row = compute_slide_family_features(slide_path=data_root / slide_name)
        feature_row["donor_id"] = donor_id
        feature_rows.append(feature_row)

    feature_df = pd.DataFrame(feature_rows)
    donor_table = cohort.merge(feature_df, on=["donor_id", "slide_name"], how="left")
    donor_table = donor_table.sort_values("donor_id").reset_index(drop=True)

    donor_feature_columns = [
        "donor_id",
        "slide_name",
        OUTCOME_COLUMN,
        "max_age_vis",
        "braak_numeric",
        "cerad_ordinal",
        "sex",
        "sex_binary",
        "mpp_um_per_px",
        "radius_px_25um",
        "n_cells_ca1",
        "n_pyramidal_ca1",
        "n_reactive_astro_ca1",
        "n_lymphocyte_ca1",
        "n_ra_exposed_pyramidal_ca1",
        "median_lymphocytes_per_exposed_pyramidal_25um",
        *(spec["feature_column"] for spec in VARIATIONS.values()),
    ]
    donor_feature_table = donor_table.loc[:, donor_feature_columns]
    donor_feature_table.to_csv("/scratch/donor_feature_table.csv", index=False)

    variation_results = [evaluate_variation(donor_table, variation_name=name) for name in VARIATIONS]
    def _rank_key(item: dict[str, Any]) -> tuple[float, float, float, float]:
        finite_sel = math.isfinite(item["selection_score"])
        finite_adj = math.isfinite(item["adjusted_score"])
        finite_loo = math.isfinite(item["loo_predictive_r"])
        return (
            1.0 if finite_sel else 0.0,
            item["selection_score"] if finite_sel else -np.inf,
            item["adjusted_score"] if finite_adj else -np.inf,
            item["loo_predictive_r"] if finite_loo else -np.inf,
        )

    ranked = sorted(variation_results, key=_rank_key, reverse=True)
    best = ranked[0]

    results_payload = {
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
        "ranked_variations": [
            {
                "variation_name": item["variation_name"],
                "description": item["description"],
                "feature_column": item["feature_column"],
                "n_total": item["n_total"],
                "n_analyzable": item["n_analyzable"],
                "coverage_fraction": item["coverage_fraction"],
                "partial_r": item["partial_r"],
                "p_value": item["p_value"],
                "selection_score": item["selection_score"],
                "loo_predictive_r": item["loo_predictive_r"],
                "is_loo_gap": item["is_loo_gap"],
                "penalty": item["penalty"],
                "adjusted_score": item["adjusted_score"],
            }
            for item in ranked
        ],
        "per_donor_loo": _json_ready(best["loo_table"].to_dict(orient="records")),
    }
    Path("/scratch/results.json").write_text(json.dumps(_json_ready(results_payload), indent=2))

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best['variation_name']}")
    print(f"  IS partial r:      {best['partial_r']:.4f}")
    print(f"  Selection score:   {best['selection_score']:.4f}")
    print(f"  LOO predictive r:  {best['loo_predictive_r']:.4f}  (diagnostic)")
    print(f"  IS-LOO Gap:        {best['is_loo_gap']:.4f}  (penalty={best['penalty']:.4f})")
    print(f"  Adjusted Score:    {best['adjusted_score']:.4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked:
        print(
            f"  {item['variation_name']}  "
            f"{item['partial_r']:.4f}  "
            f"{item['selection_score']:.4f}  "
            f"{item['loo_predictive_r']:.4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best['feature_column']}  n_ra_exposed_pyramidal_ca1")
    for _, row in best["loo_table"].iterrows():
        print(
            f"  {row['donor_id']}  "
            f"{row['outcome']:.4f}  "
            f"{row['predicted']:.4f}  "
            f"{row[best['feature_column']]:.4f}  "
            f"{int(row['n_ra_exposed_pyramidal_ca1'])}"
        )


if __name__ == "__main__":
    main()
