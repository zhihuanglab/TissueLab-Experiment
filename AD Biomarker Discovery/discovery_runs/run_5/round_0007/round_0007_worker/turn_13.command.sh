cat > /scratch/result.py <<'PY'
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import importlib.util
import json
import math
import re
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

if importlib.util.find_spec("shared_analysis") is None:  # evaluator/sandbox fallback
    shared_lib = Path("/shared/lib")
    if shared_lib.exists():
        sys.path.insert(0, str(shared_lib))

from shared_analysis import build_cell_table, load_training_cohort
from shared_analysis.artifacts import write_donor_feature_table
from shared_analysis.stats import partial_correlation


HYPOTHESIS_FAMILY = "ca1_peripyramidal_reactive_astro_lymphocyte_contact"

# These constants are rewritten by main() after the local winner is selected.
FEATURE_NAME = "__AUTO_FEATURE_NAME__"
FEATURE_COLUMN = "__AUTO_FEATURE_COLUMN__"
CANONICAL_VARIATION = "__AUTO_VARIATION__"

REGION_NAME = "CA1"
REACTIVE_ASTROCYTE_TYPE = "Reactive Astrocyte"
PYRAMIDAL_TYPE = "Pyramidal Neuron"
LYMPHOCYTE_TYPE = "Lymphocyte"

PERIPYRAMIDAL_RADIUS_UM = 30.0
DEFAULT_MPP_UM = 0.50325
MIN_REACTIVE_ASTROCYTES = 10
MIN_PYRAMIDAL_NEURONS = 25
MIN_PERIPYRAMIDAL_REACTIVE_ASTROCYTES = 10

VARIATIONS: list[dict[str, object]] = [
    {
        "name": "candidate_variant_a",
        "description": (
            "Baseline: fraction of CA1 peripyramidal reactive astrocytes with at "
            "least one CA1 lymphocyte within 30 um."
        ),
        "lymphocyte_radius_um": 30.0,
        "feature_name": "ca1_peripyramidal_reactive_astro_lymphocyte_contact_fraction_30um",
        "feature_column": "ca1_peripyramidal_reactive_astro_lymphocyte_contact_fraction_30um",
    },
    {
        "name": "candidate_variant_b",
        "description": (
            "More local version: fraction of CA1 peripyramidal reactive astrocytes "
            "with at least one CA1 lymphocyte within 20 um."
        ),
        "lymphocyte_radius_um": 20.0,
        "feature_name": "ca1_peripyramidal_reactive_astro_lymphocyte_contact_fraction_20um",
        "feature_column": "ca1_peripyramidal_reactive_astro_lymphocyte_contact_fraction_20um",
    },
]
VARIATION_BY_NAME = {spec["name"]: spec for spec in VARIATIONS}
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
OUTCOME_COLUMN = "slope_zmem0"


@lru_cache(maxsize=4)
def _load_cohort_cached(data_root_str: str) -> pd.DataFrame:
    cohort = load_training_cohort(Path(data_root_str)).copy()
    if "sex_binary" not in cohort.columns:
        sex_map = {"Female": 1.0, "Male": 0.0}
        cohort["sex_binary"] = cohort["sex"].map(sex_map)
    return cohort


def _resolve_slide_path(*, donor_row: pd.Series, data_root: Path) -> Path:
    slide_name = str(donor_row["slide_name"])
    slide_path = data_root / slide_name
    return slide_path


def _load_slide_mpp_um(slide_path: Path) -> float:
    svs_path = slide_path.with_suffix("")  # *.svs.zarr -> *.svs
    if svs_path.exists():
        try:
            import openslide  # type: ignore

            slide = openslide.OpenSlide(str(svs_path))
            for key in ("openslide.mpp-x", "aperio.MPP", "aperio.MPPX"):
                raw = slide.properties.get(key)
                if raw is None:
                    continue
                value = float(raw)
                if np.isfinite(value) and value > 0:
                    return value
        except Exception:
            pass
    return DEFAULT_MPP_UM


def _safe_tree_query(points_a: np.ndarray, points_b: np.ndarray) -> np.ndarray:
    if len(points_a) == 0 or len(points_b) == 0:
        return np.empty((0,), dtype=float)
    tree = cKDTree(points_b)
    distances, _ = tree.query(points_a, k=1)
    return np.asarray(distances, dtype=float)


def _extract_family_features_from_slide(slide_path: Path) -> dict[str, float]:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)

    ca1 = cells.loc[cells["region"] == REGION_NAME, ["x", "y", "cell_type"]]
    reactive_xy = ca1.loc[ca1["cell_type"] == REACTIVE_ASTROCYTE_TYPE, ["x", "y"]].to_numpy(
        dtype=float, copy=False
    )
    pyramidal_xy = ca1.loc[ca1["cell_type"] == PYRAMIDAL_TYPE, ["x", "y"]].to_numpy(
        dtype=float, copy=False
    )
    lymphocyte_xy = ca1.loc[ca1["cell_type"] == LYMPHOCYTE_TYPE, ["x", "y"]].to_numpy(
        dtype=float, copy=False
    )

    out: dict[str, float] = {
        "n_ca1_cells": float(len(ca1)),
        "n_ca1_reactive_astrocytes": float(len(reactive_xy)),
        "n_ca1_pyramidal_neurons": float(len(pyramidal_xy)),
        "n_ca1_lymphocytes": float(len(lymphocyte_xy)),
    }
    for spec in VARIATIONS:
        out[str(spec["feature_column"])] = float("nan")

    if len(reactive_xy) < MIN_REACTIVE_ASTROCYTES or len(pyramidal_xy) < MIN_PYRAMIDAL_NEURONS:
        out["n_peripyramidal_reactive_astrocytes"] = float("nan")
        return out

    mpp_um = _load_slide_mpp_um(slide_path)
    peripyramidal_radius_px = PERIPYRAMIDAL_RADIUS_UM / mpp_um
    reactive_to_pyramidal = _safe_tree_query(reactive_xy, pyramidal_xy)
    peripyramidal_mask = reactive_to_pyramidal <= peripyramidal_radius_px
    peripyramidal_xy = reactive_xy[peripyramidal_mask]
    out["n_peripyramidal_reactive_astrocytes"] = float(len(peripyramidal_xy))

    if len(peripyramidal_xy) < MIN_PERIPYRAMIDAL_REACTIVE_ASTROCYTES:
        return out

    if len(lymphocyte_xy) == 0:
        for spec in VARIATIONS:
            out[str(spec["feature_column"])] = 0.0
        return out

    peripyramidal_to_lymph = _safe_tree_query(peripyramidal_xy, lymphocyte_xy)
    for spec in VARIATIONS:
        lymph_radius_px = float(spec["lymphocyte_radius_um"]) / mpp_um
        value = float(np.mean(peripyramidal_to_lymph <= lymph_radius_px))
        out[str(spec["feature_column"])] = value
    return out


def _extract_family_features_for_donor(*, donor_row: pd.Series, data_root: Path) -> dict[str, float | str]:
    slide_path = _resolve_slide_path(donor_row=donor_row, data_root=data_root)
    metrics = _extract_family_features_from_slide(slide_path)
    result: dict[str, float | str] = {
        "donor_id": str(donor_row["donor_id"]),
        "slide_name": str(donor_row["slide_name"]),
        OUTCOME_COLUMN: float(donor_row[OUTCOME_COLUMN]),
        "max_age_vis": float(donor_row["max_age_vis"]),
        "braak_numeric": float(donor_row["braak_numeric"]),
        "cerad_ordinal": float(donor_row["cerad_ordinal"]),
        "sex_binary": float(donor_row["sex_binary"]),
    }
    result.update(metrics)
    return result


def _design_matrix(confound_df: pd.DataFrame) -> np.ndarray:
    values = confound_df.astype(float).to_numpy(dtype=float, copy=False)
    intercept = np.ones((len(confound_df), 1), dtype=float)
    return np.concatenate([intercept, values], axis=1)


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return float("nan")
    a2 = a[mask]
    b2 = b[mask]
    if np.std(a2) == 0 or np.std(b2) == 0:
        return float("nan")
    return float(np.corrcoef(a2, b2)[0, 1])


def _loo_prediction_table(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
    id_col: str = "donor_id",
) -> tuple[pd.DataFrame, float]:
    frame = df[[id_col, feature_col, outcome_col, *confounds]].dropna().reset_index(drop=True)
    rows: list[dict[str, float | str]] = []
    pred_resid_all: list[float] = []
    actual_resid_all: list[float] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train = _design_matrix(train[confounds])
        x_test = _design_matrix(test[confounds])

        y_feature_train = train[feature_col].to_numpy(dtype=float)
        y_outcome_train = train[outcome_col].to_numpy(dtype=float)
        y_feature_test = test[feature_col].to_numpy(dtype=float)
        y_outcome_test = test[outcome_col].to_numpy(dtype=float)

        beta_feature, *_ = np.linalg.lstsq(x_train, y_feature_train, rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train, y_outcome_train, rcond=None)

        resid_feature_train = y_feature_train - x_train @ beta_feature
        resid_outcome_train = y_outcome_train - x_train @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            predicted_resid = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            resid_feature_test = y_feature_test - x_test @ beta_feature
            predicted_resid = float(slope * resid_feature_test[0])

        actual_resid = float(y_outcome_test[0] - (x_test @ beta_outcome)[0])
        predicted_raw = float((x_test @ beta_outcome)[0] + predicted_resid) if np.isfinite(predicted_resid) else float("nan")

        rows.append(
            {
                id_col: str(test.loc[0, id_col]),
                "outcome": float(y_outcome_test[0]),
                "predicted": predicted_raw,
                "actual_resid": actual_resid,
                "predicted_resid": predicted_resid,
                feature_col: float(y_feature_test[0]),
            }
        )
        pred_resid_all.append(predicted_resid)
        actual_resid_all.append(actual_resid)

    loo_r = _corr(np.asarray(pred_resid_all, dtype=float), np.asarray(actual_resid_all, dtype=float))
    return pd.DataFrame(rows), loo_r


def _gap_penalty(is_partial_r: float, loo_predictive_r: float) -> tuple[float, float]:
    gap = float(abs(is_partial_r - loo_predictive_r)) if np.isfinite(loo_predictive_r) else float("inf")
    if not np.isfinite(gap):
        return gap, 1.0
    penalty = max(0.0, gap - 0.15) * 0.5
    return gap, penalty


def _evaluate_variation(df: pd.DataFrame, *, variation: dict[str, object], n_total: int) -> dict[str, object]:
    feature_col = str(variation["feature_column"])
    analyzable = df[["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]].dropna().copy()
    pc = partial_correlation(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    n_analyzable = int(pc["n"])
    partial_r = float(pc["partial_r"])
    selection_score = float(abs(partial_r) * (n_analyzable / n_total)) if np.isfinite(partial_r) else float("nan")
    loo_table, loo_predictive_r = _loo_prediction_table(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        id_col="donor_id",
    )
    loo_table = loo_table.merge(
        analyzable[
            [
                "donor_id",
                feature_col,
                "n_ca1_lymphocytes",
                "n_peripyramidal_reactive_astrocytes",
                "n_ca1_reactive_astrocytes",
                "n_ca1_pyramidal_neurons",
            ]
        ],
        on="donor_id",
        how="left",
        suffixes=("", "_dup"),
    )
    gap, penalty = _gap_penalty(partial_r, loo_predictive_r)
    adjusted_score = selection_score - penalty if np.isfinite(selection_score) and np.isfinite(penalty) else float("nan")
    if np.isfinite(gap) and gap > 0.30:
        adjusted_score = -1.0

    return {
        "variation_name": str(variation["name"]),
        "description": str(variation["description"]),
        "feature_name": str(variation["feature_name"]),
        "feature_column": feature_col,
        "lymphocyte_radius_um": float(variation["lymphocyte_radius_um"]),
        "n_analyzable": n_analyzable,
        "n_total": int(n_total),
        "coverage": float(n_analyzable / n_total),
        "partial_r": partial_r,
        "selection_score": selection_score,
        "loo_predictive_r": float(loo_predictive_r),
        "is_loo_gap": gap,
        "penalty": penalty,
        "adjusted_score": adjusted_score,
        "per_donor_loo": loo_table.to_dict(orient="records"),
    }


def _rank_variations(metrics: list[dict[str, object]]) -> list[dict[str, object]]:
    def key_fn(item: dict[str, object]):
        sel = float(item["selection_score"]) if np.isfinite(float(item["selection_score"])) else -np.inf
        adj = float(item["adjusted_score"]) if np.isfinite(float(item["adjusted_score"])) else -np.inf
        loo = abs(float(item["loo_predictive_r"])) if np.isfinite(float(item["loo_predictive_r"])) else -np.inf
        return (sel, adj, loo)

    return sorted(metrics, key=key_fn, reverse=True)


def _canonical_from_variation(variation_name: str) -> tuple[str, str]:
    spec = VARIATION_BY_NAME[variation_name]
    return str(spec["feature_name"]), str(spec["feature_column"])


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Canonical replay target for the family-local winner.

    Population:
      CA1 Reactive Astrocytes.
    Niche:
      subset with at least one CA1 Pyramidal Neuron within 30 um
      ("peripyramidal reactive astrocytes").
    Scalar:
      fraction of those niche-defined reactive astrocytes with at least one
      CA1 Lymphocyte within the winning radius.
    Observable pattern:
      more CA1 reactive astrocytes sitting in a neuron-adjacent inflammatory
      microenvironment.
    """
    if CANONICAL_VARIATION not in VARIATION_BY_NAME:
        raise RuntimeError(
            "Canonical variation has not been materialized yet. "
            "Run this script once as __main__ to select and rewrite the winner."
        )

    data_root = Path(data_root)
    cohort = _load_cohort_cached(str(data_root))
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    donor_row = donor_rows.iloc[0]
    features = _extract_family_features_for_donor(donor_row=donor_row, data_root=data_root)
    feature_name, feature_col = _canonical_from_variation(CANONICAL_VARIATION)
    _ = feature_name  # explicit for readability / replay auditability
    score = features.get(feature_col, float("nan"))
    if score is None or not np.isfinite(float(score)):
        return None
    return float(score)


def _rewrite_canonical_constants(best: dict[str, object]) -> None:
    path = Path(__file__)
    text = path.read_text()
    feature_name = str(best["feature_name"])
    feature_column = str(best["feature_column"])
    variation_name = str(best["variation_name"])

    replacements = {
        r'^FEATURE_NAME = ".*"$': f'FEATURE_NAME = "{feature_name}"',
        r'^FEATURE_COLUMN = ".*"$': f'FEATURE_COLUMN = "{feature_column}"',
        r'^CANONICAL_VARIATION = ".*"$': f'CANONICAL_VARIATION = "{variation_name}"',
    }
    updated = text
    for pattern, replacement in replacements.items():
        updated = re.sub(pattern, replacement, updated, flags=re.MULTILINE)
    if updated != text:
        path.write_text(updated)


def _build_donor_feature_frame(data_root: Path) -> pd.DataFrame:
    cohort = _load_cohort_cached(str(data_root))
    rows = [_extract_family_features_for_donor(donor_row=row, data_root=data_root) for _, row in cohort.iterrows()]
    return pd.DataFrame(rows)


def _print_report(best: dict[str, object], ranked: list[dict[str, object]]) -> None:
    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best['variation_name']}")
    print(f"  IS partial r:      {float(best['partial_r']):.4f}")
    print(f"  Selection score:   {float(best['selection_score']):.4f}")
    print(f"  LOO predictive r:  {float(best['loo_predictive_r']):.4f}  (diagnostic)")
    print(
        f"  IS-LOO Gap:        {float(best['is_loo_gap']):.4f}  "
        f"(penalty={float(best['penalty']):.4f})"
    )
    print(f"  Adjusted Score:    {float(best['adjusted_score']):.4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked:
        print(
            f"  {item['variation_name']}  "
            f"{float(item['partial_r']):.4f}  "
            f"{float(item['selection_score']):.4f}  "
            f"{float(item['loo_predictive_r']):.4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    best_feature = str(best["feature_column"])
    print(
        f"  donor_id  outcome  predicted  {best_feature}  "
        "n_ca1_lymphocytes  n_peripyramidal_reactive_astrocytes"
    )
    for row in best["per_donor_loo"]:
        print(
            f"  {row['donor_id']}  "
            f"{float(row['outcome']):.4f}  "
            f"{float(row['predicted']):.4f}  "
            f"{float(row[best_feature]):.6f}  "
            f"{float(row['n_ca1_lymphocytes']):.0f}  "
            f"{float(row['n_peripyramidal_reactive_astrocytes']):.0f}"
        )


def main() -> None:
    data_root = Path("/data")
    donor_df = _build_donor_feature_frame(data_root)
    n_total = int(len(donor_df))

    metrics = [_evaluate_variation(donor_df, variation=variation, n_total=n_total) for variation in VARIATIONS]
    ranked = _rank_variations(metrics)
    best = ranked[0]

    _rewrite_canonical_constants(best)

    best_feature_column = str(best["feature_column"])
    write_donor_feature_table(
        "/scratch/donor_feature_table.csv",
        donor_df,
        feature_column=best_feature_column,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUND_COLUMNS,
        id_columns=["donor_id", "slide_name"],
        extra_columns=[
            "n_ca1_cells",
            "n_ca1_reactive_astrocytes",
            "n_ca1_pyramidal_neurons",
            "n_ca1_lymphocytes",
            "n_peripyramidal_reactive_astrocytes",
        ],
    )

    payload = {
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": str(best["variation_name"]),
        "feature_name": str(best["feature_name"]),
        "feature_column": best_feature_column,
        "partial_r": float(best["partial_r"]),
        "selection_score": float(best["selection_score"]),
        "loo_predictive_r": float(best["loo_predictive_r"]),
        "is_loo_gap": float(best["is_loo_gap"]),
        "penalty": float(best["penalty"]),
        "adjusted_score": float(best["adjusted_score"]),
        "n_analyzable": int(best["n_analyzable"]),
        "n_total": int(best["n_total"]),
        "ranked_variations": [
            {
                "variation_name": str(item["variation_name"]),
                "description": str(item["description"]),
                "feature_name": str(item["feature_name"]),
                "feature_column": str(item["feature_column"]),
                "lymphocyte_radius_um": float(item["lymphocyte_radius_um"]),
                "n_analyzable": int(item["n_analyzable"]),
                "n_total": int(item["n_total"]),
                "coverage": float(item["coverage"]),
                "partial_r": float(item["partial_r"]),
                "selection_score": float(item["selection_score"]),
                "loo_predictive_r": float(item["loo_predictive_r"]),
                "is_loo_gap": float(item["is_loo_gap"]),
                "penalty": float(item["penalty"]),
                "adjusted_score": float(item["adjusted_score"]),
            }
            for item in ranked
        ],
        "per_donor_loo": best["per_donor_loo"],
        "biological_interpretation": {
            "population": "CA1 Reactive Astrocytes",
            "niche": "peripyramidal subset defined by a CA1 Pyramidal Neuron within 30 um",
            "scalar": "fraction of niche-defined reactive astrocytes with at least one CA1 Lymphocyte within the winning radius",
            "observable_pattern": "more neuron-adjacent reactive astrocytes embedded in a local inflammatory lymphocyte niche",
        },
    }
    Path("/scratch/results.json").write_text(json.dumps(payload, indent=2))
    _print_report(best, ranked)


if __name__ == "__main__":
    main()
PY
python /scratch/result.py > /scratch/run.log 2>&1
status=$?
echo "script exit status: $status"
sed -n '1,220p' /scratch/run.log
exit $status
