from __future__ import annotations

import json
import math
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, "/shared/lib")

from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort
from shared_analysis.stats import partial_correlation


HYPOTHESIS_FAMILY = "ca1_reactive_astro_fraction"
OUTCOME_COL = "slope_zmem0"
CONFOUND_COLS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
ASTROCYTE_LABEL = "Astrocyte"
REACTIVE_ASTROCYTE_LABEL = "Reactive Astrocyte"
SIDEcar_RESULTS = "results.json"

warnings.filterwarnings(
    "ignore",
    message="Object at .* is not recognized as a component of a Zarr hierarchy.",
)


@dataclass(frozen=True)
class VariationSpec:
    name: str
    regions: tuple[str, ...]
    feature_column: str
    numerator_column: str
    denominator_column: str
    region_cells_column: str
    description: str


VARIATIONS: dict[str, VariationSpec] = {
    "candidate_variant_a": VariationSpec(
        name="candidate_variant_a",
        regions=("CA1",),
        feature_column="ca1_reactive_astro_fraction",
        numerator_column="ca1_reactive_astro_n",
        denominator_column="ca1_astro_lineage_n",
        region_cells_column="ca1_region_cells_n",
        description="Reactive Astrocyte / (Astrocyte + Reactive Astrocyte) within CA1",
    ),
    "candidate_variant_b": VariationSpec(
        name="candidate_variant_b",
        regions=("CA1", "CA2"),
        feature_column="ca1_ca2_reactive_astro_fraction",
        numerator_column="ca1_ca2_reactive_astro_n",
        denominator_column="ca1_ca2_astro_lineage_n",
        region_cells_column="ca1_ca2_region_cells_n",
        description="Reactive Astrocyte / (Astrocyte + Reactive Astrocyte) within CA1+CA2",
    ),
}
BASELINE_VARIATION = "candidate_variant_a"


def _design_matrix(frame: pd.DataFrame) -> np.ndarray:
    x = frame.to_numpy(dtype=float)
    intercept = np.ones((len(frame), 1), dtype=float)
    return np.concatenate([intercept, x], axis=1)


def _compute_fraction_from_cells(cells: pd.DataFrame, regions: tuple[str, ...]) -> dict[str, float]:
    region_mask = cells["region"].isin(regions)
    region_cells = cells.loc[region_mask]
    region_cell_count = float(len(region_cells))
    lineage = region_cells.loc[
        region_cells["cell_type"].isin([ASTROCYTE_LABEL, REACTIVE_ASTROCYTE_LABEL])
    ]
    reactive_n = float((lineage["cell_type"] == REACTIVE_ASTROCYTE_LABEL).sum())
    denom_n = float(len(lineage))
    score = reactive_n / denom_n if denom_n > 0 else float("nan")
    return {
        "score": score,
        "reactive_n": reactive_n,
        "denom_n": denom_n,
        "region_cells_n": region_cell_count,
    }


def _load_selected_variation_from_sidecar() -> str:
    sidecar = Path(__file__).with_name(SIDEcar_RESULTS)
    if sidecar.exists():
        try:
            payload = json.loads(sidecar.read_text())
            name = payload.get("best_variation")
            if name in VARIATIONS:
                return str(name)
        except Exception:
            pass
    return BASELINE_VARIATION


def _extract_slide_variation_values(slide_path: Path) -> dict[str, float]:
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    out: dict[str, float] = {}
    for name, spec in VARIATIONS.items():
        comp = _compute_fraction_from_cells(cells, spec.regions)
        out[spec.feature_column] = comp["score"]
        out[spec.numerator_column] = comp["reactive_n"]
        out[spec.denominator_column] = comp["denom_n"]
        out[spec.region_cells_column] = comp["region_cells_n"]
    return out


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Canonical replay target: reactive astrocyte fraction for the selected winner.
    The selected variation is read from results.json next to this script.
    """
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    cohort = cohort.copy()
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    values = _extract_slide_variation_values(slide_path)
    variation_name = _load_selected_variation_from_sidecar()
    feature_column = VARIATIONS[variation_name].feature_column
    value = values.get(feature_column, float("nan"))
    if not np.isfinite(value):
        return None
    return float(value)


def extract_all_features(data_root: str | Path = "/data") -> pd.DataFrame:
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = (cohort["sex"].astype(str).str.lower() == "male").astype(float)

    rows: list[dict[str, Any]] = []
    for _, row in cohort.iterrows():
        slide_path = data_root / str(row["slide_name"])
        feature_values = _extract_slide_variation_values(slide_path)
        record: dict[str, Any] = {
            "donor_id": row["donor_id"],
            "slide_name": row["slide_name"],
            OUTCOME_COL: float(row[OUTCOME_COL]),
            "max_age_vis": float(row["max_age_vis"]),
            "braak_numeric": float(row["braak_numeric"]),
            "cerad_ordinal": float(row["cerad_ordinal"]),
            "sex_binary": float(row["sex_binary"]),
            "sex": row["sex"],
        }
        record.update(feature_values)
        rows.append(record)
    return pd.DataFrame(rows)


def _loo_predictions(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
    id_col: str = "donor_id",
) -> tuple[float, pd.DataFrame]:
    cols = [id_col, feature_col, outcome_col, *confounds]
    frame = df.loc[:, cols].dropna().reset_index(drop=True)
    records: list[dict[str, float | str]] = []
    pred_resids: list[float] = []
    actual_resids: list[float] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train = _design_matrix(train[confounds])
        x_test = _design_matrix(test[confounds])

        beta_feature, *_ = np.linalg.lstsq(
            x_train, train[feature_col].to_numpy(dtype=float), rcond=None
        )
        beta_outcome, *_ = np.linalg.lstsq(
            x_train, train[outcome_col].to_numpy(dtype=float), rcond=None
        )

        feature_train = train[feature_col].to_numpy(dtype=float)
        outcome_train = train[outcome_col].to_numpy(dtype=float)
        feature_test = test[feature_col].to_numpy(dtype=float)
        outcome_test = test[outcome_col].to_numpy(dtype=float)

        conf_feature_train = x_train @ beta_feature
        conf_feature_test = x_test @ beta_feature
        conf_outcome_test = x_test @ beta_outcome

        resid_feature_train = feature_train - conf_feature_train
        resid_outcome_train = outcome_train - (x_train @ beta_outcome)

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0 or not np.isfinite(denom):
            pred_resid = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            pred_resid = float(slope * (feature_test - conf_feature_test)[0])

        actual_resid = float(outcome_test[0] - conf_outcome_test[0])
        pred_raw = float(conf_outcome_test[0] + pred_resid) if np.isfinite(pred_resid) else float("nan")

        pred_resids.append(pred_resid)
        actual_resids.append(actual_resid)
        records.append(
            {
                "donor_id": str(test.loc[0, id_col]),
                "outcome": float(outcome_test[0]),
                "predicted": pred_raw,
                "predicted_residual": pred_resid,
                "actual_residual": actual_resid,
                feature_col: float(feature_test[0]),
            }
        )

    pred_arr = np.asarray(pred_resids, dtype=float)
    act_arr = np.asarray(actual_resids, dtype=float)
    mask = np.isfinite(pred_arr) & np.isfinite(act_arr)
    if mask.sum() < 3 or np.std(pred_arr[mask]) == 0 or np.std(act_arr[mask]) == 0:
        loo_r = float("nan")
    else:
        loo_r = float(np.corrcoef(pred_arr[mask], act_arr[mask])[0, 1])
    return loo_r, pd.DataFrame.from_records(records)


def _gap_penalty(gap: float) -> float:
    if not np.isfinite(gap):
        return float("nan")
    return float(max(0.0, gap - 0.15) * 0.5)


def evaluate_variation(features: pd.DataFrame, variation_name: str) -> dict[str, Any]:
    spec = VARIATIONS[variation_name]
    usable = features.dropna(subset=[spec.feature_column, OUTCOME_COL, *CONFOUND_COLS]).copy()
    n_total = int(len(features))
    n_analyzable = int(len(usable))

    metrics = partial_correlation(
        usable,
        feature_col=spec.feature_column,
        outcome_col=OUTCOME_COL,
        confounds=CONFOUND_COLS,
    )
    partial_r = float(metrics["partial_r"]) if n_analyzable >= 3 else float("nan")
    p_value = float(metrics["p_value"]) if n_analyzable >= 3 else float("nan")
    selection_score = abs(partial_r) * (n_analyzable / n_total) if np.isfinite(partial_r) else float("nan")

    loo_r, per_donor = _loo_predictions(
        usable,
        feature_col=spec.feature_column,
        outcome_col=OUTCOME_COL,
        confounds=CONFOUND_COLS,
    )
    gap = abs(abs(partial_r) - loo_r) if np.isfinite(partial_r) and np.isfinite(loo_r) else float("nan")
    penalty = _gap_penalty(gap)
    adjusted = selection_score - penalty if np.isfinite(selection_score) and np.isfinite(penalty) else float("nan")

    per_donor_base = per_donor.drop(columns=[spec.feature_column], errors="ignore")
    merged = usable.merge(per_donor_base, on="donor_id", how="left")
    residual = merged["outcome"] - merged["predicted"]
    merged["abs_error"] = residual.abs()

    return {
        "variation_name": variation_name,
        "description": spec.description,
        "regions": list(spec.regions),
        "feature_column": spec.feature_column,
        "numerator_column": spec.numerator_column,
        "denominator_column": spec.denominator_column,
        "region_cells_column": spec.region_cells_column,
        "partial_r": partial_r,
        "selection_score": float(selection_score),
        "loo_predictive_r": float(loo_r),
        "gap": float(gap),
        "gap_penalty": float(penalty),
        "adjusted_score": float(adjusted),
        "p_value": p_value,
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "per_donor": merged[
            [
                "donor_id",
                OUTCOME_COL,
                spec.feature_column,
                spec.numerator_column,
                spec.denominator_column,
                spec.region_cells_column,
                "predicted",
                "predicted_residual",
                "actual_residual",
                "abs_error",
            ]
        ]
        .rename(columns={OUTCOME_COL: "outcome"})
        .sort_values("donor_id")
        .to_dict(orient="records"),
    }


def _sort_key(result: dict[str, Any]) -> tuple[float, float, float]:
    primary = result["selection_score"]
    secondary = result["adjusted_score"]
    tertiary = result["loo_predictive_r"]
    if not np.isfinite(primary):
        primary = -np.inf
    if not np.isfinite(secondary):
        secondary = -np.inf
    if not np.isfinite(tertiary):
        tertiary = -np.inf
    return (primary, secondary, tertiary)


def run_sweep(data_root: str | Path = "/data") -> dict[str, Any]:
    features = extract_all_features(data_root)
    results = [evaluate_variation(features, name) for name in VARIATIONS]
    ranked = sorted(results, key=_sort_key, reverse=True)
    best = ranked[0]

    donor_feature_cols = [
        "donor_id",
        "slide_name",
        OUTCOME_COL,
        *CONFOUND_COLS,
        best["feature_column"],
        best["numerator_column"],
        best["denominator_column"],
        best["region_cells_column"],
    ]
    features.loc[:, donor_feature_cols].to_csv("/scratch/donor_feature_table.csv", index=False)

    payload = {
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": best["variation_name"],
        "feature_column": best["feature_column"],
        "partial_r": best["partial_r"],
        "selection_score": best["selection_score"],
        "loo_predictive_r": best["loo_predictive_r"],
        "gap": best["gap"],
        "gap_penalty": best["gap_penalty"],
        "adjusted_score": best["adjusted_score"],
        "n_total": best["n_total"],
        "n_analyzable": best["n_analyzable"],
        "ranked_variations": [
            {
                "variation_name": res["variation_name"],
                "description": res["description"],
                "regions": res["regions"],
                "feature_column": res["feature_column"],
                "partial_r": res["partial_r"],
                "selection_score": res["selection_score"],
                "loo_predictive_r": res["loo_predictive_r"],
                "gap": res["gap"],
                "gap_penalty": res["gap_penalty"],
                "adjusted_score": res["adjusted_score"],
                "n_analyzable": res["n_analyzable"],
                "n_total": res["n_total"],
                "p_value": res["p_value"],
            }
            for res in ranked
        ],
        "best_variation_spec": {
            "name": best["variation_name"],
            "description": best["description"],
            "regions": best["regions"],
            "feature_column": best["feature_column"],
            "formula": "Reactive Astrocyte / (Astrocyte + Reactive Astrocyte)",
            "population": "astrocyte-lineage cells",
        },
        "per_donor": best["per_donor"],
    }
    Path("/scratch/results.json").write_text(json.dumps(payload, indent=2))

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best['variation_name']}")
    print(f"  IS partial r:      {best['partial_r']:.4f}")
    print(f"  Selection score:   {best['selection_score']:.4f}")
    print(f"  LOO predictive r:  {best['loo_predictive_r']:.4f}  (diagnostic)")
    print(f"  IS-LOO Gap:        {best['gap']:.4f}  (penalty={best['gap_penalty']:.4f})")
    print(f"  Adjusted Score:    {best['adjusted_score']:.4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for res in ranked:
        print(
            f"  {res['variation_name']}  "
            f"{res['partial_r']:.4f}  "
            f"{res['selection_score']:.4f}  "
            f"{res['loo_predictive_r']:.4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    for row in best["per_donor"]:
        print(
            f"  {row['donor_id']}  {row['outcome']:.6f}  {row['predicted']:.6f}  "
            f"{best['feature_column']}={row[best['feature_column']]:.4f}; "
            f"{best['numerator_column']}={int(row[best['numerator_column']])}; "
            f"{best['denominator_column']}={int(row[best['denominator_column']])}; "
            f"{best['region_cells_column']}={int(row[best['region_cells_column']])}"
        )

    return {
        "features": features,
        "ranked": ranked,
        "best": best,
        "payload": payload,
    }


if __name__ == "__main__":
    run_sweep("/data")
