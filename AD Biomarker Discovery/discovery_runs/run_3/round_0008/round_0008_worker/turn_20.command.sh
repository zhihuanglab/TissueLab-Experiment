set -euo pipefail
cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.append('/shared/lib')
from shared_analysis import (  # noqa: E402
    build_cell_table,
    load_training_cohort,
    partial_correlation,
    write_donor_feature_table,
)

try:  # noqa: E402
    import openslide
except Exception:  # pragma: no cover
    openslide = None


DATA_ROOT = Path("/data")
CONTEXT_PATH = Path("/scratch/context_bundle.json")
RESULTS_PATH = Path("/scratch/results.json")
REPORT_PATH = Path("/scratch/report.md")
DONOR_FEATURE_TABLE_PATH = Path("/scratch/donor_feature_table.csv")
SELECTED_VARIATION_PATH = Path(__file__).with_name("best_variation.json")

HYPOTHESIS_FAMILY = "ca1_pyramidal_immune_admixed_reactive_cuff_fraction"
OUTCOME_COL = "slope_zmem0"
CONFOUNDS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
FALLBACK_MPP_UM = 0.503

VARIATIONS = {
    "ca1_pyramidal_immune_admixed_reactive_cuff_fraction_r35um_ra2_ly1": {
        "radius_um": 35.0,
        "reactive_min": 2,
        "lymph_min": 1,
        "feature_column": "ca1_pyramidal_immune_admixed_reactive_cuff_fraction_r35um_ra2_ly1",
        "description": "Fraction of CA1 pyramidal neurons with >=2 Reactive Astrocytes and >=1 Lymphocyte within 35 um.",
    },
    "ca1_pyramidal_immune_admixed_reactive_cuff_fraction_r50um_ra2_ly1": {
        "radius_um": 50.0,
        "reactive_min": 2,
        "lymph_min": 1,
        "feature_column": "ca1_pyramidal_immune_admixed_reactive_cuff_fraction_r50um_ra2_ly1",
        "description": "Fraction of CA1 pyramidal neurons with >=2 Reactive Astrocytes and >=1 Lymphocyte within 50 um.",
    },
}
BASELINE_VARIATION = "ca1_pyramidal_immune_admixed_reactive_cuff_fraction_r35um_ra2_ly1"
DEFAULT_CANONICAL_VARIATION = BASELINE_VARIATION


def _load_worker_brief() -> dict:
    if not CONTEXT_PATH.exists():
        return {}
    obj = json.loads(CONTEXT_PATH.read_text())
    return obj.get("worker_brief", {})


def _selected_variation_name() -> str:
    if SELECTED_VARIATION_PATH.exists():
        try:
            payload = json.loads(SELECTED_VARIATION_PATH.read_text())
            name = payload.get("best_variation")
            if name in VARIATIONS:
                return name
        except Exception:
            pass
    return DEFAULT_CANONICAL_VARIATION


FEATURE_NAME = _selected_variation_name()
FEATURE_COLUMN = VARIATIONS[FEATURE_NAME]["feature_column"]


def _safe_float(value) -> float | None:
    if value is None:
        return None
    try:
        value = float(value)
    except Exception:
        return None
    if not math.isfinite(value):
        return None
    return value


def _json_ready(value):
    if isinstance(value, dict):
        return {str(k): _json_ready(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_json_ready(v) for v in value]
    if isinstance(value, tuple):
        return [_json_ready(v) for v in value]
    if isinstance(value, (np.floating, np.integer)):
        value = value.item()
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return None
    return value


def _design_matrix(frame: pd.DataFrame) -> np.ndarray:
    cols = [np.ones(len(frame), dtype=float)]
    for col in frame.columns:
        cols.append(frame[col].to_numpy(dtype=float))
    return np.column_stack(cols)


def _get_slide_mpp_um(slide_stem: Path) -> float:
    svs_path = slide_stem.with_suffix("")
    if not str(svs_path).endswith(".svs"):
        svs_path = Path(str(slide_stem).replace(".svs.zarr", ".svs"))
    if openslide is None or not svs_path.exists():
        return FALLBACK_MPP_UM
    try:
        slide = openslide.OpenSlide(str(svs_path))
        props = slide.properties
        for key in ("openslide.mpp-x", "aperio.MPP", "aperio.MPPX"):
            raw = props.get(key)
            if raw is not None:
                value = float(raw)
                if math.isfinite(value) and value > 0:
                    slide.close()
                    return value
        slide.close()
    except Exception:
        return FALLBACK_MPP_UM
    return FALLBACK_MPP_UM


def _compute_variation_scores_for_cells(
    ca1_cells: pd.DataFrame,
    *,
    mpp_um: float,
) -> dict[str, float]:
    pyramidal = ca1_cells.loc[ca1_cells["cell_type"] == "Pyramidal Neuron", ["x", "y"]].to_numpy(dtype=float)
    reactive = ca1_cells.loc[ca1_cells["cell_type"] == "Reactive Astrocyte", ["x", "y"]].to_numpy(dtype=float)
    lymph = ca1_cells.loc[ca1_cells["cell_type"] == "Lymphocyte", ["x", "y"]].to_numpy(dtype=float)

    if len(pyramidal) == 0:
        return {name: float("nan") for name in VARIATIONS}

    reactive_tree = cKDTree(reactive) if len(reactive) else None
    lymph_tree = cKDTree(lymph) if len(lymph) else None

    scores: dict[str, float] = {}
    for name, spec in VARIATIONS.items():
        radius_px = float(spec["radius_um"]) / float(mpp_um)
        reactive_counts = (
            reactive_tree.query_ball_point(pyramidal, r=radius_px, return_length=True)
            if reactive_tree is not None
            else np.zeros(len(pyramidal), dtype=int)
        )
        lymph_counts = (
            lymph_tree.query_ball_point(pyramidal, r=radius_px, return_length=True)
            if lymph_tree is not None
            else np.zeros(len(pyramidal), dtype=int)
        )
        cuffed = (np.asarray(reactive_counts) >= int(spec["reactive_min"])) & (
            np.asarray(lymph_counts) >= int(spec["lymph_min"])
        )
        scores[name] = float(np.mean(cuffed))
    return scores


def _extract_donor_row(donor_row: pd.Series, *, data_root: Path) -> dict:
    slide_name = str(donor_row["slide_name"])
    slide_path = data_root / slide_name
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)

    ca1_cells = cells.loc[cells["region"] == "CA1", ["x", "y", "cell_type"]].copy()
    mpp_um = _get_slide_mpp_um(slide_path)
    variation_scores = _compute_variation_scores_for_cells(ca1_cells, mpp_um=mpp_um)

    pyramidal_n = int((ca1_cells["cell_type"] == "Pyramidal Neuron").sum())
    reactive_n = int((ca1_cells["cell_type"] == "Reactive Astrocyte").sum())
    lymph_n = int((ca1_cells["cell_type"] == "Lymphocyte").sum())

    row = donor_row.to_dict()
    row["sex_binary"] = 1.0 if str(donor_row.get("sex", "")).strip().lower() == "male" else 0.0
    row["ca1_pyramidal_n"] = pyramidal_n
    row["ca1_reactive_astrocyte_n"] = reactive_n
    row["ca1_lymphocyte_n"] = lymph_n
    row["mpp_um"] = float(mpp_um)
    for name, value in variation_scores.items():
        row[VARIATIONS[name]["feature_column"]] = value
    return row


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    row = _extract_donor_row(donor_rows.iloc[0], data_root=data_root)
    selected_name = _selected_variation_name()
    return row.get(VARIATIONS[selected_name]["feature_column"])


def _loo_prediction_table(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
) -> tuple[pd.DataFrame, float]:
    frame = df.loc[:, ["donor_id", feature_col, outcome_col, *confounds]].dropna().reset_index(drop=True)
    if len(frame) < 3:
        return frame.assign(predicted=np.nan, predicted_residual=np.nan, outcome_residual=np.nan), float("nan")

    preds_raw: list[float] = []
    preds_resid: list[float] = []
    acts_resid: list[float] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train_conf = _design_matrix(train[confounds])
        x_test_conf = _design_matrix(test[confounds])

        y_feature_train = train[feature_col].to_numpy(dtype=float)
        y_outcome_train = train[outcome_col].to_numpy(dtype=float)
        y_outcome_test = test[outcome_col].to_numpy(dtype=float)

        beta_feature, *_ = np.linalg.lstsq(x_train_conf, y_feature_train, rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train_conf, y_outcome_train, rcond=None)

        resid_feature_train = y_feature_train - x_train_conf @ beta_feature
        resid_outcome_train = y_outcome_train - x_train_conf @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            preds_raw.append(float("nan"))
            preds_resid.append(float("nan"))
            acts_resid.append(float("nan"))
            continue

        slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
        resid_feature_test = test[feature_col].to_numpy(dtype=float) - x_test_conf @ beta_feature
        outcome_resid_test = y_outcome_test - x_test_conf @ beta_outcome
        pred_resid = float(slope * resid_feature_test[0])
        pred_raw = float((x_test_conf @ beta_outcome)[0] + pred_resid)

        preds_raw.append(pred_raw)
        preds_resid.append(pred_resid)
        acts_resid.append(float(outcome_resid_test[0]))

    out = frame.copy()
    out["predicted"] = preds_raw
    out["predicted_residual"] = preds_resid
    out["outcome_residual"] = acts_resid

    mask = np.isfinite(out["predicted_residual"].to_numpy()) & np.isfinite(out["outcome_residual"].to_numpy())
    if mask.sum() >= 3 and np.std(out.loc[mask, "predicted_residual"]) > 0 and np.std(out.loc[mask, "outcome_residual"]) > 0:
        loo_r = float(np.corrcoef(out.loc[mask, "predicted_residual"], out.loc[mask, "outcome_residual"])[0, 1])
    else:
        loo_r = float("nan")
    return out, loo_r


def _evaluate_variation(df: pd.DataFrame, variation_name: str) -> dict:
    feature_col = VARIATIONS[variation_name]["feature_column"]
    analyzable = df.loc[:, [feature_col, OUTCOME_COL, *CONFOUNDS]].dropna()
    pcorr = partial_correlation(
        df,
        feature_col=feature_col,
        outcome_col=OUTCOME_COL,
        confounds=CONFOUNDS,
    )
    pred_table, loo_r = _loo_prediction_table(
        df,
        feature_col=feature_col,
        outcome_col=OUTCOME_COL,
        confounds=CONFOUNDS,
    )
    n_total = int(len(df))
    n_analyzable = int(len(analyzable))
    partial_r = float(pcorr["partial_r"])
    selection_score = abs(partial_r) * (n_analyzable / n_total) if n_total else float("nan")
    gap = abs(partial_r) - loo_r if math.isfinite(partial_r) and math.isfinite(loo_r) else float("nan")
    penalty = gap if math.isfinite(gap) else float("nan")
    adjusted_score = selection_score - penalty if math.isfinite(selection_score) and math.isfinite(penalty) else float("nan")
    return {
        "variation_name": variation_name,
        "feature_column": feature_col,
        "description": VARIATIONS[variation_name]["description"],
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "partial_r": partial_r,
        "selection_score": float(selection_score),
        "loo_predictive_r": float(loo_r),
        "is_loo_gap": float(gap),
        "penalty": float(penalty),
        "adjusted_score": float(adjusted_score),
        "p_value": float(pcorr["p_value"]) if math.isfinite(pcorr["p_value"]) else float("nan"),
        "predictions": pred_table,
    }


def _format_float(value: float | None) -> str:
    value = _safe_float(value)
    return "NA" if value is None else f"{value:.4f}"


def main() -> None:
    brief = _load_worker_brief()
    cohort = load_training_cohort(DATA_ROOT).copy()
    donor_rows = []
    for _, donor_row in cohort.iterrows():
        donor_rows.append(_extract_donor_row(donor_row, data_root=DATA_ROOT))
    donor_df = pd.DataFrame(donor_rows)

    results = [_evaluate_variation(donor_df, name) for name in VARIATIONS]
    ranked = sorted(
        results,
        key=lambda r: (
            -np.nan_to_num(r["selection_score"], nan=-np.inf),
            -np.nan_to_num(r["adjusted_score"], nan=-np.inf),
        ),
    )
    best = ranked[0]

    SELECTED_VARIATION_PATH.write_text(
        json.dumps(
            {
                "best_variation": best["variation_name"],
                "feature_column": best["feature_column"],
                "hypothesis_family": HYPOTHESIS_FAMILY,
            },
            indent=2,
        )
    )

    extra_columns = [
        col
        for col in [
            "ca1_pyramidal_n",
            "ca1_reactive_astrocyte_n",
            "ca1_lymphocyte_n",
            "mpp_um",
            *(VARIATIONS[name]["feature_column"] for name in VARIATIONS if VARIATIONS[name]["feature_column"] != best["feature_column"]),
        ]
        if col in donor_df.columns
    ]
    write_donor_feature_table(
        DONOR_FEATURE_TABLE_PATH,
        donor_df,
        feature_column=best["feature_column"],
        outcome_column=OUTCOME_COL,
        covariates=CONFOUNDS,
        extra_columns=extra_columns,
    )

    ranked_payload = []
    for item in ranked:
        ranked_payload.append(
            {
                k: _json_ready(v)
                for k, v in item.items()
                if k != "predictions"
            }
        )

    predictions_table = best["predictions"].copy()
    predictions_table = predictions_table.merge(
        donor_df[["donor_id", best["feature_column"], "ca1_pyramidal_n", "ca1_reactive_astrocyte_n", "ca1_lymphocyte_n"]],
        on="donor_id",
        how="left",
    )

    results_payload = {
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": best["variation_name"],
        "feature_name": best["variation_name"],
        "feature_column": best["feature_column"],
        "n_total": best["n_total"],
        "n_analyzable": best["n_analyzable"],
        "partial_r": _json_ready(best["partial_r"]),
        "selection_score": _json_ready(best["selection_score"]),
        "loo_predictive_r": _json_ready(best["loo_predictive_r"]),
        "is_loo_gap": _json_ready(best["is_loo_gap"]),
        "penalty": _json_ready(best["penalty"]),
        "adjusted_score": _json_ready(best["adjusted_score"]),
        "p_value": _json_ready(best["p_value"]),
        "ranked_variations": ranked_payload,
        "artifacts": {
            "donor_feature_table": str(DONOR_FEATURE_TABLE_PATH),
            "selected_variation": str(SELECTED_VARIATION_PATH),
        },
        "per_donor_loo": _json_ready(
            predictions_table[
                [
                    "donor_id",
                    OUTCOME_COL,
                    "predicted",
                    "predicted_residual",
                    "outcome_residual",
                    best["feature_column"],
                    "ca1_pyramidal_n",
                    "ca1_reactive_astrocyte_n",
                    "ca1_lymphocyte_n",
                ]
            ].to_dict(orient="records")
        ),
    }
    RESULTS_PATH.write_text(json.dumps(_json_ready(results_payload), indent=2))

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
    for item in ranked:
        print(
            f"  {item['variation_name']}  "
            f"{_format_float(item['partial_r'])}  "
            f"{_format_float(item['selection_score'])}  "
            f"{_format_float(item['loo_predictive_r'])}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best['feature_column']}")
    for _, row in predictions_table.sort_values("donor_id").iterrows():
        print(
            f"  {row['donor_id']}  "
            f"{_format_float(row[OUTCOME_COL])}  "
            f"{_format_float(row['predicted'])}  "
            f"{_format_float(row[best['feature_column']])}"
        )


if __name__ == "__main__":
    main()
PY
python /scratch/result.py > /scratch/run_output.txt
sed -n '1,220p' /scratch/run_output.txt
