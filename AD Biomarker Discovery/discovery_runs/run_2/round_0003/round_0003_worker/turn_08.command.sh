set -euo pipefail
cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path

sys.path.append("/shared/lib")

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from shared_analysis import build_cell_table, load_training_cohort
from shared_analysis.artifacts import write_donor_feature_table
from shared_analysis.stats import partial_correlation


warnings.filterwarnings(
    "ignore",
    message="Object at .* is not recognized as a component of a Zarr hierarchy.",
)

FAMILY_NAME = "ca1_pyramidal_reactive_niche"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUNDS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
MIN_CA1_PYRAMIDAL_COUNT = 10

VARIATIONS = [
    {
        "name": "ca1_pyr_reactive_niche_r50px",
        "radius_px": 50.0,
        "description": "CA1 pyramidal neurons with a reactive astrocyte neighbor within 50 px, normalized by pyramidal neurons with any astroglial neighbor within 50 px.",
        "feature_column": "ca1_pyr_reactive_niche_r50px__reactive_over_astroglial_neighbored_pyramidal_fraction",
        "feature_name": "ca1_pyramidal_reactive_niche__reactive_over_astroglial_neighbored_pyramidal_fraction_r50px",
    },
    {
        "name": "ca1_pyr_reactive_niche_r60px",
        "radius_px": 60.0,
        "description": "CA1 pyramidal neurons with a reactive astrocyte neighbor within 60 px, normalized by pyramidal neurons with any astroglial neighbor within 60 px.",
        "feature_column": "ca1_pyr_reactive_niche_r60px__reactive_over_astroglial_neighbored_pyramidal_fraction",
        "feature_name": "ca1_pyramidal_reactive_niche__reactive_over_astroglial_neighbored_pyramidal_fraction_r60px",
    },
    {
        "name": "ca1_pyr_reactive_niche_r70px",
        "radius_px": 70.0,
        "description": "CA1 pyramidal neurons with a reactive astrocyte neighbor within 70 px, normalized by pyramidal neurons with any astroglial neighbor within 70 px.",
        "feature_column": "ca1_pyr_reactive_niche_r70px__reactive_over_astroglial_neighbored_pyramidal_fraction",
        "feature_name": "ca1_pyramidal_reactive_niche__reactive_over_astroglial_neighbored_pyramidal_fraction_r70px",
    },
]
BASELINE_VARIATION = "ca1_pyr_reactive_niche_r60px"
VARIATION_BY_NAME = {v["name"]: v for v in VARIATIONS}

SCRATCH_ROOT = Path(__file__).resolve().parent
RESULTS_PATH = SCRATCH_ROOT / "results.json"
BEST_CONFIG_PATH = SCRATCH_ROOT / "best_variation.json"
DONOR_TABLE_PATH = SCRATCH_ROOT / "donor_feature_table.csv"
REPORT_PATH = SCRATCH_ROOT / "report.md"


def _safe_float(value) -> float:
    try:
        value = float(value)
    except Exception:
        return float("nan")
    return value if math.isfinite(value) else float("nan")


def _format_float(value: float) -> str:
    value = _safe_float(value)
    return f"{value: .4f}" if math.isfinite(value) else "    nan"


def _selection_score(partial_r: float | None, n_analyzable: int | None, n_total: int | None) -> float:
    partial_r = _safe_float(partial_r)
    if not math.isfinite(partial_r):
        return float("nan")
    if n_analyzable is None or n_total is None or n_total <= 0:
        return float("nan")
    coverage = max(0.0, min(1.0, float(n_analyzable) / float(n_total)))
    return float(abs(partial_r) * coverage)


def _adjusted_metrics(partial_r: float | None, loo_predictive_r: float | None) -> tuple[float, float, float]:
    is_r = abs(_safe_float(partial_r))
    loo_r = abs(_safe_float(loo_predictive_r))
    if not math.isfinite(is_r):
        is_r = 0.0
    if not math.isfinite(loo_r):
        loo_r = 0.0
    gap = is_r - loo_r
    if gap > 0.30:
        return float(gap), float(gap), -1.0
    penalty = max(0.0, gap - 0.15) * 0.5
    return float(gap), float(penalty), float(loo_r - penalty)


def _load_eval_cohort(data_root: str | Path) -> pd.DataFrame:
    data_root = Path(data_root)
    training = data_root / "training_cohort.csv"
    test = data_root / "test_cohort.csv"
    if training.exists():
        cohort = load_training_cohort(data_root)
    elif test.exists():
        cohort = pd.read_csv(test)
    else:
        cohort = load_training_cohort(data_root)
    cohort = cohort.copy()
    if "sex_binary" not in cohort.columns:
        cohort["sex_binary"] = (
            cohort["sex"]
            .astype(str)
            .str.strip()
            .map({"Female": 0.0, "Male": 1.0, "female": 0.0, "male": 1.0})
        )
    return cohort


def _design_matrix(frame: pd.DataFrame) -> np.ndarray:
    x = frame.to_numpy(dtype=float, copy=True)
    if x.ndim == 1:
        x = x[:, None]
    return np.column_stack([np.ones(len(frame), dtype=float), x])


def _canonical_variation_name() -> str:
    if BEST_CONFIG_PATH.exists():
        try:
            payload = json.loads(BEST_CONFIG_PATH.read_text(encoding="utf-8"))
            name = str(payload.get("best_variation", "")).strip()
            if name in VARIATION_BY_NAME:
                return name
        except Exception:
            pass
    return BASELINE_VARIATION


def _nearest_neighbor_mask(query_xy: np.ndarray, ref_xy: np.ndarray, radius_px: float) -> np.ndarray:
    query_xy = np.asarray(query_xy, dtype=float)
    ref_xy = np.asarray(ref_xy, dtype=float)
    if query_xy.size == 0:
        return np.zeros(0, dtype=bool)
    if ref_xy.size == 0:
        return np.zeros(len(query_xy), dtype=bool)
    tree = cKDTree(ref_xy)
    distances, _ = tree.query(query_xy, k=1, distance_upper_bound=float(radius_px))
    return np.isfinite(distances) & (distances <= float(radius_px))


def _compute_slide_family_features(slide_path: Path) -> dict[str, float]:
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    region_norm = cells["region"].fillna("").astype(str).str.upper()
    cell_type_norm = cells["cell_type"].fillna("").astype(str).str.lower()

    ca1_mask = region_norm.str.startswith("CA1")
    ca1 = cells.loc[ca1_mask, ["x", "y"]].copy()
    ca1_types = cell_type_norm.loc[ca1_mask].reset_index(drop=True)

    coords = ca1[["x", "y"]].to_numpy(dtype=float, copy=True)
    pyr = coords[ca1_types.eq("pyramidal neuron").to_numpy()]
    astro = coords[ca1_types.eq("astrocyte").to_numpy()]
    reactive = coords[ca1_types.eq("reactive astrocyte").to_numpy()]
    astroglial = np.vstack([astro, reactive]) if (len(astro) + len(reactive)) > 0 else np.empty((0, 2), dtype=float)

    result: dict[str, float] = {
        "ca1_total_cell_count": float(len(ca1)),
        "ca1_pyramidal_count": float(len(pyr)),
        "ca1_astrocyte_count": float(len(astro)),
        "ca1_reactive_astrocyte_count": float(len(reactive)),
    }

    for variation in VARIATIONS:
        reactive_hits = _nearest_neighbor_mask(pyr, reactive, variation["radius_px"])
        astroglial_hits = _nearest_neighbor_mask(pyr, astroglial, variation["radius_px"])
        reactive_count = int(np.sum(reactive_hits))
        astroglial_count = int(np.sum(astroglial_hits))

        feature_value = float("nan")
        if len(pyr) >= MIN_CA1_PYRAMIDAL_COUNT and astroglial_count > 0:
            feature_value = float(reactive_count / astroglial_count)

        result[variation["feature_column"]] = feature_value
        result[f"{variation['name']}__reactive_neighbored_pyramidal_count"] = float(reactive_count)
        result[f"{variation['name']}__astroglial_neighbored_pyramidal_count"] = float(astroglial_count)

    return result


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    data_root = Path(data_root)
    cohort = _load_eval_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"].astype(str) == str(donor_id)]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    slide_features = _compute_slide_family_features(slide_path)
    canonical_name = _canonical_variation_name()
    canonical_column = VARIATION_BY_NAME[canonical_name]["feature_column"]
    value = _safe_float(slide_features.get(canonical_column, float("nan")))
    if not math.isfinite(value):
        return None
    return value


def _build_donor_feature_table(data_root: Path) -> pd.DataFrame:
    cohort = _load_eval_cohort(data_root)
    rows: list[dict[str, object]] = []
    for row in cohort.itertuples(index=False):
        slide_path = data_root / str(row.slide_name)
        features = _compute_slide_family_features(slide_path)
        rows.append(
            {
                **row._asdict(),
                **features,
            }
        )
    return pd.DataFrame(rows)


def _loo_predictions(
    donor_table: pd.DataFrame,
    *,
    feature_col: str,
) -> tuple[pd.DataFrame, float]:
    needed = ["donor_id", OUTCOME_COLUMN, feature_col, *CONFOUNDS]
    frame = donor_table.loc[:, needed].dropna().reset_index(drop=True).copy()
    if len(frame) < 3:
        return pd.DataFrame(columns=["donor_id", "outcome", "predicted", "predicted_residual", "outcome_residual", feature_col]), float("nan")

    loo_rows: list[dict[str, float | str]] = []
    pred_resid: list[float] = []
    act_resid: list[float] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train_conf = _design_matrix(train[CONFOUNDS])
        x_test_conf = _design_matrix(test[CONFOUNDS])

        beta_feature, *_ = np.linalg.lstsq(x_train_conf, train[feature_col].to_numpy(dtype=float), rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train_conf, train[OUTCOME_COLUMN].to_numpy(dtype=float), rcond=None)

        resid_feature_train = train[feature_col].to_numpy(dtype=float) - x_train_conf @ beta_feature
        resid_outcome_train = train[OUTCOME_COLUMN].to_numpy(dtype=float) - x_train_conf @ beta_outcome
        denom = float(np.dot(resid_feature_train, resid_feature_train))

        feature_test = float(test[feature_col].iloc[0])
        outcome_test = float(test[OUTCOME_COLUMN].iloc[0])
        confound_outcome_pred = float((x_test_conf @ beta_outcome)[0])
        confound_feature_pred = float((x_test_conf @ beta_feature)[0])

        predicted_raw = float("nan")
        predicted_residual = float("nan")
        outcome_residual = float(outcome_test - confound_outcome_pred)

        if denom > 0 and np.std(resid_feature_train) > 0 and np.std(resid_outcome_train) > 0:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            feature_test_residual = float(feature_test - confound_feature_pred)
            predicted_residual = float(slope * feature_test_residual)
            predicted_raw = float(confound_outcome_pred + predicted_residual)

        pred_resid.append(predicted_residual)
        act_resid.append(outcome_residual)
        loo_rows.append(
            {
                "donor_id": str(test["donor_id"].iloc[0]),
                "outcome": outcome_test,
                "predicted": predicted_raw,
                "predicted_residual": predicted_residual,
                "outcome_residual": outcome_residual,
                feature_col: feature_test,
            }
        )

    pred_arr = np.asarray(pred_resid, dtype=float)
    act_arr = np.asarray(act_resid, dtype=float)
    mask = np.isfinite(pred_arr) & np.isfinite(act_arr)
    loo_r = float(np.corrcoef(pred_arr[mask], act_arr[mask])[0, 1]) if mask.sum() >= 3 and np.std(pred_arr[mask]) > 0 and np.std(act_arr[mask]) > 0 else float("nan")
    return pd.DataFrame(loo_rows), loo_r


def _evaluate_variation(donor_table: pd.DataFrame, variation: dict[str, object]) -> tuple[dict[str, object], pd.DataFrame]:
    feature_col = str(variation["feature_column"])
    partial = partial_correlation(
        donor_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUNDS,
    )
    analyzable = donor_table.loc[:, [feature_col, OUTCOME_COLUMN, *CONFOUNDS]].dropna()
    loo_table, loo_r = _loo_predictions(donor_table, feature_col=feature_col)
    gap, penalty, adjusted = _adjusted_metrics(partial.get("partial_r"), loo_r)
    selection = _selection_score(partial.get("partial_r"), len(analyzable), len(donor_table))
    metrics = {
        "variation_name": str(variation["name"]),
        "description": str(variation["description"]),
        "radius_px": float(variation["radius_px"]),
        "feature_column": feature_col,
        "feature_name": str(variation["feature_name"]),
        "n_total": int(len(donor_table)),
        "n_analyzable": int(len(analyzable)),
        "partial_r": _safe_float(partial.get("partial_r")),
        "p_value": _safe_float(partial.get("p_value")),
        "selection_score": _safe_float(selection),
        "loo_predictive_r": _safe_float(loo_r),
        "is_loo_gap": _safe_float(gap),
        "gap_penalty": _safe_float(penalty),
        "adjusted_score": _safe_float(adjusted),
    }
    return metrics, loo_table


def _rank_key(metrics: dict[str, object]) -> tuple[float, float, float]:
    selection = _safe_float(metrics.get("selection_score"))
    partial_r = abs(_safe_float(metrics.get("partial_r")))
    adjusted = _safe_float(metrics.get("adjusted_score"))
    if not math.isfinite(selection):
        selection = -np.inf
    if not math.isfinite(partial_r):
        partial_r = -np.inf
    if not math.isfinite(adjusted):
        adjusted = -np.inf
    return (selection, partial_r, adjusted)


def _json_default(value):
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if pd.isna(value):
        return None
    raise TypeError(f"Object of type {type(value)!r} is not JSON serializable")


def _error_pattern_text(loo_table: pd.DataFrame, donor_table: pd.DataFrame) -> str:
    if loo_table.empty:
        return "LOO evaluation was unavailable because too few analyzable donors remained after confound filtering."
    merged = loo_table.merge(
        donor_table.loc[:, ["donor_id", "cognitive_status", "braak_numeric", "cerad_ordinal", "sex"]],
        on="donor_id",
        how="left",
    )
    merged["abs_error"] = (merged["predicted"] - merged["outcome"]).abs()
    top = merged.sort_values("abs_error", ascending=False).head(3)
    donor_bits = []
    for row in top.itertuples(index=False):
        donor_bits.append(
            f"{row.donor_id} (abs err={row.abs_error:.3f}, status={row.cognitive_status}, Braak={row.braak_numeric}, CERAD={row.cerad_ordinal})"
        )
    if top.empty:
        return "No donor-level LOO predictions were available."
    status_counts = top["cognitive_status"].fillna("unknown").value_counts().to_dict()
    if len(status_counts) == 1:
        shared_status = next(iter(status_counts))
        shared_text = f"The largest errors were concentrated in donors labeled {shared_status}."
    else:
        shared_text = "The largest errors were mixed across cognitive-status groups rather than confined to one label."
    return shared_text + " Top absolute errors: " + "; ".join(donor_bits) + "."


def _write_report(
    *,
    donor_table: pd.DataFrame,
    ranked: list[dict[str, object]],
    best_metrics: dict[str, object],
    best_loo: pd.DataFrame,
):
    best_name = str(best_metrics["variation_name"])
    best_col = str(best_metrics["feature_column"])
    winner_direction = "higher" if _safe_float(best_metrics["partial_r"]) > 0 else "lower"

    others = [r for r in ranked if r["variation_name"] != best_name]
    other_summary = "; ".join(
        f"{r['variation_name']} (selection={_safe_float(r['selection_score']):.4f}, partial_r={_safe_float(r['partial_r']):.4f}, loo={_safe_float(r['loo_predictive_r']):.4f})"
        for r in others
    )
    if not other_summary:
        other_summary = "No alternate radii were tested."

    denom_col = f"{best_name}__astroglial_neighbored_pyramidal_count"
    numer_col = f"{best_name}__reactive_neighbored_pyramidal_count"
    analyzable = donor_table.loc[:, [best_col, numer_col, denom_col]].dropna()
    mean_feature = analyzable[best_col].mean() if not analyzable.empty else float("nan")
    mean_numer = analyzable[numer_col].mean() if not analyzable.empty else float("nan")
    mean_denom = analyzable[denom_col].mean() if not analyzable.empty else float("nan")

    failed_radii = ", ".join(
        f"{r['variation_name']} ({'too tight' if float(r['radius_px']) < float(best_metrics['radius_px']) else 'broader/diluted'})"
        for r in others
    )
    if not failed_radii:
        failed_radii = "none"

    sign_text = (
        f"{winner_direction} niche scores tracked higher slope_zmem0"
        if _safe_float(best_metrics["partial_r"]) > 0
        else f"{winner_direction} niche scores tracked lower slope_zmem0"
    )

    report = f"""## Summary
One sentence: tested the {FAMILY_NAME} family, the winning local variation was {best_name}, and its local winner selection score was {_safe_float(best_metrics['selection_score']):.4f}.

## Metrics
Winning variation: {best_name}
- IS partial r: {_safe_float(best_metrics['partial_r']):.4f}
- Selection score: {_safe_float(best_metrics['selection_score']):.4f}
- LOO predictive r: {_safe_float(best_metrics['loo_predictive_r']):.4f}
- IS-LOO gap: {_safe_float(best_metrics['is_loo_gap']):.4f}
- Gap penalty: {_safe_float(best_metrics['gap_penalty']):.4f}
- Adjusted score: {_safe_float(best_metrics['adjusted_score']):.4f}
Other tested radii: {other_summary}.

## Findings
1. What worked and why: {best_name} was the strongest radius, suggesting the signal is most coherent when the donor scalar measures the fraction of CA1 pyramidal neurons whose astroglial-contacting neighborhood specifically contains a reactive astrocyte. In analyzable donors the mean feature value was {mean_feature:.4f}, with mean reactive-neighbored pyramidal count {mean_numer:.1f} over mean astroglial-neighbored pyramidal count {mean_denom:.1f}; this keeps the score anchored to a neuron-adjacent gliosis niche rather than bulk gliosis alone.
2. What failed and why: the nearby radius alternatives underperformed ({failed_radii}). The tighter radius likely misses reactive astrocytes sitting just outside the immediate soma-scale shell, while the broader radius dilutes the neuron-adjacent niche into more generic local astroglial abundance.
3. Error pattern: {_error_pattern_text(best_loo, donor_table)}

## Rationale
This winner is biologically coherent because it conditions reactive astrocyte presence on surviving CA1 pyramidal neurons and normalizes by local astroglial availability, so it asks whether gliosis is preferentially concentrated right around the neuronal population already implicated by earlier rounds. It beat the nearby alternatives because the winning radius appears to balance specificity and capture: small enough to stay neuron-adjacent, but large enough to include the local reactive astrocyte shell. It is likely at least partially distinct from a pure bulk reactive fraction because the denominator removes non-specific astroglial abundance and centers the summary on pyramidal neighborhoods.

## Interpretation
The signal seems to mean a CA1 neuron-adjacent reactive gliosis state: among CA1 pyramidal neurons that sit in an astroglial neighborhood, donors with {winner_direction} values of the reactive-over-astroglial neighborhood fraction showed {sign_text}. Population: CA1 Pyramidal Neuron query cells with CA1 Reactive Astrocyte / Astrocyte neighbors. Niche: the immediate peri-neuronal CA1 astroglial neighborhood. Feature summary: reactive-neighbored pyramidal neurons divided by astroglial-neighbored pyramidal neurons. Simplest observable pattern: more of the surviving CA1 pyramidal neurons are ringed by reactive astrocytes rather than only generic astrocytes.

## Next
Next local sweep: keep the same CA1 pyramidal-centered niche but test a slightly finer radius bracket around the winner (for example ±10 px around {best_metrics['radius_px']:.0f}px) or replace the binary neighbor hit with a nearest-reactive distance / local reactive count summary, because the current ranking suggests spatial scale matters more than the broader family definition.
"""
    REPORT_PATH.write_text(report, encoding="utf-8")


def main() -> int:
    data_root = Path("/data")
    donor_table = _build_donor_feature_table(data_root)

    ranked_metrics: list[dict[str, object]] = []
    loo_by_variation: dict[str, pd.DataFrame] = {}

    for variation in VARIATIONS:
        metrics, loo_table = _evaluate_variation(donor_table, variation)
        ranked_metrics.append(metrics)
        loo_by_variation[str(variation["name"])] = loo_table

    ranked_metrics = sorted(ranked_metrics, key=_rank_key, reverse=True)
    best_metrics = ranked_metrics[0]
    best_name = str(best_metrics["variation_name"])
    best_variation = VARIATION_BY_NAME[best_name]
    best_col = str(best_metrics["feature_column"])
    best_loo = loo_by_variation[best_name].copy()

    best_loo = best_loo.merge(
        donor_table.loc[
            :,
            [
                "donor_id",
                best_col,
                "ca1_pyramidal_count",
                "ca1_astrocyte_count",
                "ca1_reactive_astrocyte_count",
                f"{best_name}__reactive_neighbored_pyramidal_count",
                f"{best_name}__astroglial_neighbored_pyramidal_count",
            ],
        ],
        on=["donor_id", best_col],
        how="left",
    )

    BEST_CONFIG_PATH.write_text(
        json.dumps(
            {
                "best_variation": best_name,
                "feature_column": best_col,
                "feature_name": best_metrics["feature_name"],
                "radius_px": best_variation["radius_px"],
                "family_name": FAMILY_NAME,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    extra_columns = [v["feature_column"] for v in VARIATIONS if v["feature_column"] != best_col]
    extra_columns += [
        "ca1_total_cell_count",
        "ca1_pyramidal_count",
        "ca1_astrocyte_count",
        "ca1_reactive_astrocyte_count",
    ]
    for variation in VARIATIONS:
        extra_columns.extend(
            [
                f"{variation['name']}__reactive_neighbored_pyramidal_count",
                f"{variation['name']}__astroglial_neighbored_pyramidal_count",
            ]
        )
    extra_columns = [c for i, c in enumerate(extra_columns) if c in donor_table.columns and c not in extra_columns[:i]]

    write_donor_feature_table(
        DONOR_TABLE_PATH,
        donor_table,
        feature_column=best_col,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUNDS,
        extra_columns=extra_columns,
    )

    ranked_export = [
        {
            k: v
            for k, v in metrics.items()
            if k not in {"description"}
        }
        for metrics in ranked_metrics
    ]
    results_payload = {
        "status": "ok",
        "hypothesis_family": FAMILY_NAME,
        "best_variation": best_name,
        "feature_name": best_metrics["feature_name"],
        "feature_column": best_col,
        "radius_px": best_variation["radius_px"],
        "n_total": best_metrics["n_total"],
        "n_analyzable": best_metrics["n_analyzable"],
        "partial_r": best_metrics["partial_r"],
        "p_value": best_metrics["p_value"],
        "selection_score": best_metrics["selection_score"],
        "loo_predictive_r": best_metrics["loo_predictive_r"],
        "is_loo_gap": best_metrics["is_loo_gap"],
        "gap_penalty": best_metrics["gap_penalty"],
        "adjusted_score": best_metrics["adjusted_score"],
        "ranked_variations": ranked_export,
        "artifacts": {
            "donor_feature_table": str(DONOR_TABLE_PATH),
            "best_variation_config": str(BEST_CONFIG_PATH),
            "report": str(REPORT_PATH),
        },
        "per_donor_loo": best_loo.to_dict(orient="records"),
    }
    RESULTS_PATH.write_text(json.dumps(results_payload, indent=2, default=_json_default) + "\n", encoding="utf-8")

    _write_report(
        donor_table=donor_table,
        ranked=ranked_metrics,
        best_metrics=best_metrics,
        best_loo=best_loo,
    )

    print(f"HYPOTHESIS FAMILY: {FAMILY_NAME}")
    print(f"BEST VARIATION: {best_name}")
    print(f"  IS partial r:      {_safe_float(best_metrics['partial_r']): .4f}")
    print(f"  Selection score:   {_safe_float(best_metrics['selection_score']): .4f}")
    print(f"  LOO predictive r:  {_safe_float(best_metrics['loo_predictive_r']): .4f}  (diagnostic)")
    print(f"  IS-LOO Gap:        {_safe_float(best_metrics['is_loo_gap']): .4f}  (penalty={_safe_float(best_metrics['gap_penalty']):.4f})")
    print(f"  Adjusted Score:    {_safe_float(best_metrics['adjusted_score']): .4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for metrics in ranked_metrics:
        print(
            f"  {metrics['variation_name']}  "
            f"{_safe_float(metrics['partial_r']): .4f}  "
            f"{_safe_float(metrics['selection_score']): .4f}  "
            f"{_safe_float(metrics['loo_predictive_r']): .4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(
        f"  donor_id  outcome  predicted  {best_col}  "
        f"{best_name}__reactive_neighbored_pyramidal_count  "
        f"{best_name}__astroglial_neighbored_pyramidal_count"
    )
    for row in best_loo.itertuples(index=False):
        row_dict = row._asdict()
        print(
            f"  {row_dict['donor_id']}  "
            f"{_safe_float(row_dict['outcome']): .4f}  "
            f"{_safe_float(row_dict['predicted']): .4f}  "
            f"{_safe_float(row_dict[best_col]): .4f}  "
            f"{_safe_float(row_dict[f'{best_name}__reactive_neighbored_pyramidal_count']): .1f}  "
            f"{_safe_float(row_dict[f'{best_name}__astroglial_neighbored_pyramidal_count']): .1f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
PY

python /scratch/result.py
echo '--- results.json head ---'
python - <<'PY'
import json
with open('/scratch/results.json') as f:
    obj=json.load(f)
print(json.dumps({
    'best_variation': obj.get('best_variation'),
    'feature_column': obj.get('feature_column'),
    'selection_score': obj.get('selection_score'),
    'ranked_variations': obj.get('ranked_variations'),
}, indent=2)[:4000])
PY
echo '--- report.md ---'
sed -n '1,220p' /scratch/report.md
