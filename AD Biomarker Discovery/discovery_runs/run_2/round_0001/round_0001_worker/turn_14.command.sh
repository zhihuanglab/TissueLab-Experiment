cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath

SHARED_LIB = Path("/shared/lib")
if SHARED_LIB.exists() and str(SHARED_LIB) not in sys.path:
    sys.path.insert(0, str(SHARED_LIB))

from shared_analysis import (  # type: ignore
    load_centroids,
    load_class_ids,
    load_class_lookup,
    load_region_polygons,
    load_training_cohort,
    partial_correlation,
    write_donor_feature_table,
)

warnings.filterwarnings("ignore", message="Object at .* is not recognized as a component of a Zarr hierarchy.")

DATA_ROOT_DEFAULT = Path("/data")
SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_PATH = SCRIPT_DIR / "results.json"
REPORT_PATH = SCRIPT_DIR / "report.md"
DONOR_TABLE_PATH = SCRIPT_DIR / "donor_feature_table.csv"
REPLAY_CONFIG_PATH = SCRIPT_DIR / "result_replay.json"

HYPOTHESIS_FAMILY = "CA1 pyramidal neuron burden"
OUTCOME_COL = "slope_zmem0"
CONFOUND_COLS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
BASELINE_VARIATION = "ca1_pyramidal_fraction"
VARIATIONS = [
    {
        "name": "ca1_pyramidal_fraction",
        "description": "CA1 Pyramidal Neuron count divided by all classified CA1 cells.",
    },
    {
        "name": "ca1_pyramidal_density",
        "description": "CA1 Pyramidal Neuron count divided by annotated CA1 area in mm^2.",
    },
]


def _load_replay_config() -> dict[str, Any]:
    if REPLAY_CONFIG_PATH.exists():
        try:
            return json.loads(REPLAY_CONFIG_PATH.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}


CANONICAL_VARIATION = str(_load_replay_config().get("best_variation", BASELINE_VARIATION))
FEATURE_NAME = f"ca1_pyramidal_burden__{CANONICAL_VARIATION}"
FEATURE_COLUMN = CANONICAL_VARIATION


def _safe_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        value = float(value)
    except Exception:
        return None
    if not math.isfinite(value):
        return None
    return value


def _fmt(value: Any) -> str:
    value = _safe_float(value)
    return "nan" if value is None else f"{value:.4f}"


def _jsonify(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _jsonify(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_jsonify(v) for v in value]
    if isinstance(value, tuple):
        return [_jsonify(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        value = float(value)
        if math.isnan(value) or math.isinf(value):
            return None
        return value
    if isinstance(value, (pd.Timestamp,)):
        return value.isoformat()
    return value


def polygon_area_px2(points: np.ndarray) -> float:
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[0] < 3:
        return float("nan")
    x = pts[:, 0]
    y = pts[:, 1]
    return float(abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))) * 0.5)


def get_slide_mpp_um(slide_path: Path) -> float:
    import openslide

    slide = openslide.OpenSlide(str(slide_path))
    prop_keys = ["openslide.mpp-x", "aperio.MPP", "openslide.mpp-y"]
    for key in prop_keys:
        raw = slide.properties.get(key)
        if raw is None:
            continue
        try:
            value = float(raw)
        except Exception:
            continue
        if math.isfinite(value) and value > 0:
            return value
    raise RuntimeError(f"Could not read MPP from {slide_path}")


def compute_variations_for_slide(*, zarr_path: Path, svs_path: Path) -> dict[str, float]:
    centroids = load_centroids(zarr_path)
    class_ids = load_class_ids(zarr_path)
    class_lookup = load_class_lookup(zarr_path)

    region_polygons = load_region_polygons(zarr_path, scale=16.0).get("CA1", [])
    if not region_polygons:
        return {
            "ca1_pyramidal_fraction": float("nan"),
            "ca1_pyramidal_density": float("nan"),
            "ca1_pyramidal_count": float("nan"),
            "ca1_total_classified_cells": float("nan"),
            "ca1_area_mm2": float("nan"),
        }

    ca1_mask = np.zeros(len(centroids), dtype=bool)
    for polygon in region_polygons:
        if len(polygon) < 3:
            continue
        ca1_mask |= MplPath(np.asarray(polygon, dtype=float), closed=True).contains_points(centroids)

    ca1_total = int(ca1_mask.sum())
    if ca1_total <= 0:
        return {
            "ca1_pyramidal_fraction": float("nan"),
            "ca1_pyramidal_density": float("nan"),
            "ca1_pyramidal_count": 0.0,
            "ca1_total_classified_cells": 0.0,
            "ca1_area_mm2": float("nan"),
        }

    pyramidal_ids = [int(cid) for cid, name in class_lookup.items() if str(name) == "Pyramidal Neuron"]
    pyramidal_mask = np.isin(class_ids, pyramidal_ids) if pyramidal_ids else np.zeros(len(class_ids), dtype=bool)
    pyramidal_count = int(np.sum(ca1_mask & pyramidal_mask))

    area_px2 = float(sum(polygon_area_px2(poly) for poly in region_polygons if len(poly) >= 3))
    mpp_um = get_slide_mpp_um(svs_path)
    area_mm2 = float(area_px2 * (mpp_um ** 2) / 1_000_000.0) if area_px2 > 0 else float("nan")

    fraction = float(pyramidal_count / ca1_total) if ca1_total > 0 else float("nan")
    density = float(pyramidal_count / area_mm2) if area_mm2 and math.isfinite(area_mm2) and area_mm2 > 0 else float("nan")

    return {
        "ca1_pyramidal_fraction": fraction,
        "ca1_pyramidal_density": density,
        "ca1_pyramidal_count": float(pyramidal_count),
        "ca1_total_classified_cells": float(ca1_total),
        "ca1_area_mm2": area_mm2,
    }


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"].astype(str) == str(donor_id)]
    if donor_rows.empty:
        return None

    slide_name = str(donor_rows.iloc[0]["slide_name"])
    zarr_path = data_root / slide_name
    svs_path = data_root / slide_name.replace(".zarr", "")
    features = compute_variations_for_slide(zarr_path=zarr_path, svs_path=svs_path)

    replay = _load_replay_config()
    canonical_variation = str(replay.get("best_variation", BASELINE_VARIATION))
    value = _safe_float(features.get(canonical_variation))
    return None if value is None else float(value)


def prepare_donor_feature_table(data_root: Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map({"Female": 0.0, "Male": 1.0})

    records: list[dict[str, Any]] = []
    for _, row in cohort.iterrows():
        slide_name = str(row["slide_name"])
        zarr_path = data_root / slide_name
        svs_path = data_root / slide_name.replace(".zarr", "")
        features = compute_variations_for_slide(zarr_path=zarr_path, svs_path=svs_path)

        record = row.to_dict()
        record.update(features)
        records.append(record)

    return pd.DataFrame.from_records(records)


def loo_prediction_table(df: pd.DataFrame, *, feature_col: str) -> tuple[pd.DataFrame, float]:
    needed_cols = ["donor_id", "slide_name", OUTCOME_COL, feature_col, *CONFOUND_COLS]
    frame = df.loc[:, needed_cols].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)

    preds_full: list[float] = []
    preds_resid: list[float] = []
    actuals_full: list[float] = []
    actuals_resid: list[float] = []
    donor_ids: list[str] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train = np.column_stack([np.ones(len(train)), train[CONFOUND_COLS].astype(float).to_numpy()])
        x_test = np.column_stack([np.ones(len(test)), test[CONFOUND_COLS].astype(float).to_numpy()])

        y_feature_train = train[feature_col].astype(float).to_numpy()
        y_outcome_train = train[OUTCOME_COL].astype(float).to_numpy()

        beta_feature, *_ = np.linalg.lstsq(x_train, y_feature_train, rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train, y_outcome_train, rcond=None)

        resid_feature_train = y_feature_train - x_train @ beta_feature
        resid_outcome_train = y_outcome_train - x_train @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            pred_resid = float("nan")
            pred_full = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            feature_test = test[feature_col].astype(float).to_numpy()
            outcome_test = test[OUTCOME_COL].astype(float).to_numpy()
            resid_feature_test = feature_test - x_test @ beta_feature
            pred_resid = float(slope * resid_feature_test[0])
            pred_full = float((x_test @ beta_outcome)[0] + pred_resid)

        outcome_test = test[OUTCOME_COL].astype(float).to_numpy()
        resid_outcome_test = outcome_test - x_test @ beta_outcome

        preds_full.append(pred_full)
        preds_resid.append(pred_resid)
        actuals_full.append(float(outcome_test[0]))
        actuals_resid.append(float(resid_outcome_test[0]))
        donor_ids.append(str(test.iloc[0]["donor_id"]))

    out = frame.copy()
    out["predicted"] = preds_full
    out["predicted_residual"] = preds_resid
    out["actual_residual"] = actuals_resid
    out["abs_error"] = np.abs(out["predicted"] - out[OUTCOME_COL])

    mask = np.isfinite(np.asarray(preds_resid)) & np.isfinite(np.asarray(actuals_resid))
    if mask.sum() >= 3 and np.std(np.asarray(preds_resid)[mask]) > 0 and np.std(np.asarray(actuals_resid)[mask]) > 0:
        loo_r = float(np.corrcoef(np.asarray(preds_resid)[mask], np.asarray(actuals_resid)[mask])[0, 1])
    else:
        loo_r = float("nan")
    return out, loo_r


def evaluate_variation(df: pd.DataFrame, variation: dict[str, str]) -> dict[str, Any]:
    feature_col = variation["name"]
    analyzable = df.loc[:, ["donor_id", OUTCOME_COL, feature_col, *CONFOUND_COLS]].replace([np.inf, -np.inf], np.nan).dropna()

    if len(analyzable) < 3:
        return {
            "variation_name": feature_col,
            "description": variation["description"],
            "feature_column": feature_col,
            "n_total": int(len(df)),
            "n_analyzable": int(len(analyzable)),
            "partial_r": float("nan"),
            "p_value": float("nan"),
            "selection_score": float("nan"),
            "loo_predictive_r": float("nan"),
            "is_loo_gap": float("nan"),
            "penalty": float("nan"),
            "adjusted_score": float("nan"),
            "loo_table": pd.DataFrame(),
        }

    pc = partial_correlation(
        df,
        feature_col=feature_col,
        outcome_col=OUTCOME_COL,
        confounds=CONFOUND_COLS,
    )
    partial_r = float(pc["partial_r"])
    p_value = float(pc["p_value"])
    n_total = int(len(df))
    n_analyzable = int(pc["n"])
    selection_score = float(abs(partial_r) * (n_analyzable / n_total)) if math.isfinite(partial_r) else float("nan")

    loo_table, loo_predictive_r = loo_prediction_table(df, feature_col=feature_col)
    is_loo_gap = float(abs(partial_r - loo_predictive_r)) if math.isfinite(partial_r) and math.isfinite(loo_predictive_r) else float("nan")
    penalty = float(max(0.0, is_loo_gap - 0.10)) if math.isfinite(is_loo_gap) else float("nan")
    adjusted_score = float(selection_score - penalty) if math.isfinite(selection_score) and math.isfinite(penalty) else float("nan")

    return {
        "variation_name": feature_col,
        "description": variation["description"],
        "feature_column": feature_col,
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "partial_r": partial_r,
        "p_value": p_value,
        "selection_score": selection_score,
        "loo_predictive_r": loo_predictive_r,
        "is_loo_gap": is_loo_gap,
        "penalty": penalty,
        "adjusted_score": adjusted_score,
        "loo_table": loo_table,
    }


def choose_best_variation(metrics_list: list[dict[str, Any]]) -> dict[str, Any]:
    order = {item["name"]: idx for idx, item in enumerate(VARIATIONS)}

    def sort_key(item: dict[str, Any]) -> tuple[float, float, int]:
        selection = _safe_float(item.get("selection_score"))
        partial_r = _safe_float(item.get("partial_r"))
        return (
            -1e9 if selection is None else selection,
            -1e9 if partial_r is None else abs(partial_r),
            -order.get(item["variation_name"], 999),
        )

    ranked = sorted(metrics_list, key=sort_key, reverse=True)
    return ranked[0]


def build_report(
    *,
    donor_table: pd.DataFrame,
    ranked_metrics: list[dict[str, Any]],
    best_metric: dict[str, Any],
) -> str:
    best_name = best_metric["variation_name"]
    winner_table = best_metric["loo_table"].copy()
    winner_feature_desc = (
        "share of all CA1 classified cells labeled Pyramidal Neuron"
        if best_name == "ca1_pyramidal_fraction"
        else "count of CA1 Pyramidal Neuron cells per mm^2 of annotated CA1"
    )
    tissue_pattern = (
        "a thinner CA1 pyramidal band relative to the total local CA1 cell population"
        if best_name == "ca1_pyramidal_fraction"
        else "sparser packing of CA1 pyramidal neurons within the annotated CA1 ribbon"
    )

    ranking_lines = []
    for idx, metric in enumerate(ranked_metrics, start=1):
        ranking_lines.append(
            f"{idx}. `{metric['variation_name']}`: partial r {_fmt(metric['partial_r'])}, "
            f"selection score {_fmt(metric['selection_score'])}, "
            f"LOO predictive r {_fmt(metric['loo_predictive_r'])}."
        )

    error_text = "No analyzable donors in LOO."
    if not winner_table.empty:
        merge_cols = ["donor_id", "cognitive_status", "braak_label", "cerad_label", "overall_ad_neuropath_change", "sex"]
        merged = winner_table.merge(donor_table[merge_cols], on="donor_id", how="left")
        worst = merged.sort_values("abs_error", ascending=False).head(5).copy()
        donor_list = ", ".join(
            f"{row['donor_id']} (|err|={row['abs_error']:.3f})"
            for _, row in worst.iterrows()
        )
        shared_bits = []
        for col, label in [
            ("cognitive_status", "cognitive status"),
            ("overall_ad_neuropath_change", "AD neuropath change"),
            ("braak_label", "Braak"),
            ("cerad_label", "CERAD"),
            ("sex", "sex"),
        ]:
            if col not in worst.columns:
                continue
            mode = worst[col].dropna().mode()
            if mode.empty:
                continue
            top_value = str(mode.iloc[0])
            hits = int((worst[col] == top_value).sum())
            if hits >= 3:
                shared_bits.append(f"{label} `{top_value}` ({hits}/5)")
        shared_text = "; ".join(shared_bits) if shared_bits else "no single clinicopathologic category dominated the worst-error donors"
        error_text = f"Worst absolute LOO errors were {donor_list}; among these donors, {shared_text}."

    if best_name == "ca1_pyramidal_fraction":
        worked_text = (
            "The fraction variant worked best because it normalizes CA1 pyramidal burden against the total classified "
            "CA1 cellular compartment, reducing slide-to-slide differences in annotation size and overall tissue amount."
        )
        failed_text = (
            "The density variant underperformed because absolute counts per annotated area appear more sensitive to "
            "regional outline size, section thickness, and local packing differences that are not specific to selective CA1 neuron loss."
        )
        rationale_text = (
            "A compositional loss signal is biologically coherent in CA1: as pyramidal neurons are depleted or replaced "
            "by glial and other non-pyramidal cells, the pyramidal share of CA1 should fall. This beat raw density because "
            "the nearby alternative retained extra geometric variance from the annotation footprint."
        )
        next_text = (
            "Next, keep the CA1 pyramidal-burden direction but test denominators that are closer to the CA1 circuit itself, "
            "for example CA1 pyramidal cells divided by all CA1 neuronal cells rather than all classified CA1 cells."
        )
    else:
        worked_text = (
            "The density variant worked best because it captures how many CA1 pyramidal neurons remain per unit annotated CA1 tissue, "
            "which is closer to an absolute depletion signal than a relative composition score."
        )
        failed_text = (
            "The fraction variant underperformed because a donor can retain a similar pyramidal share even when overall CA1 cellularity "
            "and the absolute number of pyramidal neurons are both reduced."
        )
        rationale_text = (
            "Absolute CA1 pyramidal packing is biologically coherent for a memory-circuit injury marker: selective neuronal depletion should "
            "present as fewer pyramidal-classified cells across the CA1 ribbon per unit tissue area. This beat the fraction alternative because "
            "the stronger signal came from absolute local scarcity, not merely compositional replacement."
        )
        next_text = (
            "Next, stay near this winner and test whether CA1 pyramidal density sharpens further when restricted to the central pyramidal band "
            "or contrasted against CA3 pyramidal density to reduce whole-slide technical variation."
        )

    report = f"""## Summary
Tested the {HYPOTHESIS_FAMILY} family with two local normalization variants; `{best_name}` won with selection score {_fmt(best_metric['selection_score'])}.

## Metrics
Winning variation: `{best_name}`.
- IS partial r: {_fmt(best_metric['partial_r'])}
- Selection score: {_fmt(best_metric['selection_score'])}
- LOO predictive r: {_fmt(best_metric['loo_predictive_r'])}
- IS-LOO gap: {_fmt(best_metric['is_loo_gap'])}
- Penalty: {_fmt(best_metric['penalty'])}
- Adjusted score: {_fmt(best_metric['adjusted_score'])}
- Analyzable donors: {best_metric['n_analyzable']}/{best_metric['n_total']}

Ranking of tested variations:
{chr(10).join(ranking_lines)}

## Findings
1. What worked and why: {worked_text}
2. What failed and why: {failed_text}
3. Error pattern: {error_text}

## Rationale
{rationale_text}
Because this is round 1 with no accepted panel, the winner is a plausible seed biomarker rather than evidence of additivity beyond an existing panel.

## Interpretation
Biologically, the signal looks like CA1 pyramidal neuron preservation versus depletion within hippocampal memory circuitry.
- Population: `Pyramidal Neuron`
- Niche/region: `CA1`
- Donor-level scalar: {winner_feature_desc}
- Simplest observable tissue pattern: {tissue_pattern}

## Next
{next_text}
"""
    return report


def main() -> None:
    global CANONICAL_VARIATION, FEATURE_NAME, FEATURE_COLUMN

    data_root = DATA_ROOT_DEFAULT
    donor_table = prepare_donor_feature_table(data_root)

    metrics_list = [evaluate_variation(donor_table, variation) for variation in VARIATIONS]
    order = {item["name"]: idx for idx, item in enumerate(VARIATIONS)}
    ranked_metrics = sorted(
        metrics_list,
        key=lambda item: (
            -1e9 if _safe_float(item.get("selection_score")) is None else _safe_float(item["selection_score"]),
            -1e9 if _safe_float(item.get("partial_r")) is None else abs(_safe_float(item["partial_r"])),
            -order.get(item["variation_name"], 999),
        ),
        reverse=True,
    )
    best_metric = ranked_metrics[0]
    best_variation = str(best_metric["variation_name"])

    CANONICAL_VARIATION = best_variation
    FEATURE_NAME = f"ca1_pyramidal_burden__{best_variation}"
    FEATURE_COLUMN = best_variation

    REPLAY_CONFIG_PATH.write_text(
        json.dumps(
            {
                "best_variation": best_variation,
                "feature_column": best_variation,
                "feature_name": FEATURE_NAME,
                "hypothesis_family": HYPOTHESIS_FAMILY,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    extra_columns = [
        col
        for col in [
            "ca1_pyramidal_fraction",
            "ca1_pyramidal_density",
            "ca1_pyramidal_count",
            "ca1_total_classified_cells",
            "ca1_area_mm2",
            "cognitive_status",
            "braak_label",
            "cerad_label",
            "overall_ad_neuropath_change",
            "sex",
        ]
        if col != best_variation and col in donor_table.columns
    ]
    write_donor_feature_table(
        DONOR_TABLE_PATH,
        donor_table,
        feature_column=best_variation,
        outcome_column=OUTCOME_COL,
        covariates=CONFOUND_COLS,
        id_columns=["donor_id", "slide_name"],
        extra_columns=extra_columns,
    )

    best_payload = {
        k: v
        for k, v in best_metric.items()
        if k != "loo_table"
    }
    results_payload = {
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": best_variation,
        "feature_name": FEATURE_NAME,
        "feature_column": best_variation,
        **best_payload,
        "ranked_variations": [
            {
                k: v
                for k, v in metric.items()
                if k != "loo_table"
            }
            for metric in ranked_metrics
        ],
    }
    RESULTS_PATH.write_text(json.dumps(_jsonify(results_payload), indent=2) + "\n", encoding="utf-8")

    report_text = build_report(
        donor_table=donor_table,
        ranked_metrics=ranked_metrics,
        best_metric=best_metric,
    )
    REPORT_PATH.write_text(report_text, encoding="utf-8")

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best_variation}")
    print(f"  IS partial r:      {_fmt(best_metric['partial_r'])}")
    print(f"  Selection score:   {_fmt(best_metric['selection_score'])}")
    print(f"  LOO predictive r:  {_fmt(best_metric['loo_predictive_r'])}  (diagnostic)")
    print(f"  IS-LOO Gap:        {_fmt(best_metric['is_loo_gap'])}  (penalty={_fmt(best_metric['penalty'])})")
    print(f"  Adjusted Score:    {_fmt(best_metric['adjusted_score'])}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for metric in ranked_metrics:
        print(
            f"  {metric['variation_name']}  "
            f"{_fmt(metric['partial_r'])}  "
            f"{_fmt(metric['selection_score'])}  "
            f"{_fmt(metric['loo_predictive_r'])}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  ca1_pyramidal_fraction  ca1_pyramidal_density")
    winner_table = best_metric["loo_table"].merge(
        donor_table[["donor_id", "ca1_pyramidal_fraction", "ca1_pyramidal_density"]],
        on="donor_id",
        how="left",
    )
    for _, row in winner_table.iterrows():
        print(
            f"  {row['donor_id']}  "
            f"{_fmt(row[OUTCOME_COL])}  "
            f"{_fmt(row['predicted'])}  "
            f"{_fmt(row['ca1_pyramidal_fraction'])}  "
            f"{_fmt(row['ca1_pyramidal_density'])}"
        )


if __name__ == "__main__":
    main()
PY
python -m py_compile /scratch/result.py
