cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import sys
import warnings
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, "/shared/lib")
from shared_analysis.artifacts import write_donor_feature_table
from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort
from shared_analysis.stats import partial_correlation


# Hypothesis family from the worker brief.
HYPOTHESIS_FAMILY = "ca1_reactive_neuron_exclusion"
FEATURE_NAME = HYPOTHESIS_FAMILY

# Fallback canonical replay target. main() writes canonical_variation.json with the
# empirically winning variation; compute_donor_score loads that sidecar first.
BEST_VARIATION = "reactive_pyramidal_desert_t10"
FEATURE_COLUMN = "ca1_reactive_pyramidal_desert_t10"

OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]

TARGET_REGION = "CA1"
REACTIVE_LABEL = "Reactive Astrocyte"
PYRAMIDAL_LABEL = "Pyramidal Neuron"
LOCAL_RADIUS_UM = 35.0
DEFAULT_MPP_UM = 0.5

VARIATIONS = {
    "reactive_pyramidal_desert_t10": {
        "threshold": 0.10,
        "feature_column": "ca1_reactive_pyramidal_desert_t10",
        "description": "Fraction of CA1 reactive astrocytes whose local 35 µm CA1 neighborhood contains <10% pyramidal neurons.",
    },
    "reactive_pyramidal_desert_t20": {
        "threshold": 0.20,
        "feature_column": "ca1_reactive_pyramidal_desert_t20",
        "description": "Fraction of CA1 reactive astrocytes whose local 35 µm CA1 neighborhood contains <20% pyramidal neurons.",
    },
}


def _suppress_zarr_sidecar_warnings() -> None:
    warnings.filterwarnings(
        "ignore",
        message=r"Object at .* is not recognized as a component of a Zarr hierarchy\.",
        category=UserWarning,
    )


def _load_eval_cohort(data_root: str | Path) -> pd.DataFrame:
    data_root = Path(data_root)
    training_path = data_root / "training_cohort.csv"
    test_path = data_root / "test_cohort.csv"
    if training_path.exists():
        return load_training_cohort(data_root)
    if test_path.exists():
        return pd.read_csv(test_path)
    return load_training_cohort(data_root)


def _safe_float(value) -> float | None:
    try:
        value = float(value)
    except Exception:
        return None
    if not math.isfinite(value):
        return None
    return value


def _json_default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        value = float(obj)
        return None if not math.isfinite(value) else value
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Not JSON serializable: {type(obj)!r}")


def _selection_score(partial_r: float | None, n_analyzable: int, n_total: int) -> float:
    if partial_r is None or not math.isfinite(partial_r) or n_total <= 0:
        return float("nan")
    return float(abs(partial_r) * (float(n_analyzable) / float(n_total)))


def _adjusted_metrics(partial_r: float | None, loo_predictive_r: float | None) -> tuple[float, float, float]:
    is_r = abs(float(partial_r)) if partial_r is not None and math.isfinite(float(partial_r)) else 0.0
    loo_r = abs(float(loo_predictive_r)) if loo_predictive_r is not None and math.isfinite(float(loo_predictive_r)) else 0.0
    gap = is_r - loo_r
    if gap > 0.30:
        return float(gap), float(gap), -1.0
    penalty = max(0.0, gap - 0.15) * 0.5
    return float(gap), float(penalty), float(loo_r - penalty)


def _design_matrix(df: pd.DataFrame) -> np.ndarray:
    return np.column_stack([np.ones(len(df), dtype=float)] + [df[col].to_numpy(dtype=float) for col in df.columns])


def _load_slide_mpp_um(svs_path: Path) -> float:
    try:
        import openslide  # type: ignore

        with openslide.OpenSlide(str(svs_path)) as slide:
            for key in ("openslide.mpp-x", "aperio.MPP"):
                if key in slide.properties:
                    value = _safe_float(slide.properties.get(key))
                    if value is not None and value > 0:
                        return float(value)
    except Exception:
        pass
    return DEFAULT_MPP_UM


def _radius_px_from_slide(svs_path: Path) -> float:
    return float(LOCAL_RADIUS_UM / _load_slide_mpp_um(svs_path))


def _local_pyramidal_desert_features(
    cells: pd.DataFrame,
    *,
    radius_px: float,
) -> tuple[dict[str, float], dict[str, int | float]]:
    ca1 = cells.loc[cells["region"] == TARGET_REGION, ["x", "y", "cell_type"]].copy()
    if ca1.empty:
        features = {spec["feature_column"]: float("nan") for spec in VARIATIONS.values()}
        return features, {
            "ca1_cell_count": 0,
            "ca1_reactive_count": 0,
            "ca1_pyramidal_count": 0,
            "valid_center_count": 0,
            "radius_px": radius_px,
        }

    coords = ca1[["x", "y"]].to_numpy(dtype=np.float32)
    cell_types = ca1["cell_type"].astype(str).to_numpy()
    is_reactive = cell_types == REACTIVE_LABEL
    is_pyramidal = cell_types == PYRAMIDAL_LABEL

    reactive_indices = np.flatnonzero(is_reactive)
    reactive_count = int(reactive_indices.size)
    features = {spec["feature_column"]: float("nan") for spec in VARIATIONS.values()}
    metadata = {
        "ca1_cell_count": int(len(ca1)),
        "ca1_reactive_count": reactive_count,
        "ca1_pyramidal_count": int(is_pyramidal.sum()),
        "valid_center_count": 0,
        "radius_px": radius_px,
    }
    if reactive_count == 0:
        return features, metadata

    tree = cKDTree(coords)
    neighbor_lists = tree.query_ball_point(coords[reactive_indices], r=float(radius_px))

    local_pyramidal_fracs: list[float] = []
    for center_idx, neighbors in zip(reactive_indices, neighbor_lists):
        neighbor_idx = np.asarray(neighbors, dtype=int)
        if neighbor_idx.size == 0:
            continue
        neighbor_idx = neighbor_idx[neighbor_idx != center_idx]
        if neighbor_idx.size == 0:
            continue
        local_pyramidal_fracs.append(float(is_pyramidal[neighbor_idx].mean()))

    metadata["valid_center_count"] = int(len(local_pyramidal_fracs))
    if not local_pyramidal_fracs:
        return features, metadata

    frac_array = np.asarray(local_pyramidal_fracs, dtype=float)
    for spec in VARIATIONS.values():
        threshold = float(spec["threshold"])
        features[spec["feature_column"]] = float(np.mean(frac_array < threshold))
    return features, metadata


def extract_donor_features(*, slide_path: str | Path, svs_path: str | Path) -> dict[str, float | int]:
    _suppress_zarr_sidecar_warnings()
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    radius_px = _radius_px_from_slide(Path(svs_path))
    features, metadata = _local_pyramidal_desert_features(cells, radius_px=radius_px)
    return {**features, **metadata}


def _load_canonical_config() -> dict[str, object]:
    path = Path(__file__).with_name("canonical_variation.json")
    if path.exists():
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            pass
    return {
        "best_variation": BEST_VARIATION,
        "feature_column": FEATURE_COLUMN,
        "hypothesis_family": HYPOTHESIS_FAMILY,
    }


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Return the canonical winning biomarker score for one donor.

    This replay target recomputes the best local variation from raw slide data.
    """
    data_root = Path(data_root)
    cohort = _load_eval_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"].astype(str) == str(donor_id)]
    if donor_rows.empty:
        return None

    canonical = _load_canonical_config()
    feature_column = str(canonical.get("feature_column", FEATURE_COLUMN))
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    svs_name = slide_name.replace(".svs.zarr", ".svs")
    svs_path = data_root / svs_name
    features = extract_donor_features(slide_path=slide_path, svs_path=svs_path)
    value = features.get(feature_column)
    return None if value is None or not math.isfinite(float(value)) else float(value)


@dataclass
class VariationEvaluation:
    name: str
    feature_column: str
    partial_r: float | None
    p_value: float | None
    n_analyzable: int
    n_total: int
    selection_score: float | None
    loo_predictive_r: float | None
    is_loo_gap: float | None
    gap_penalty: float | None
    adjusted_score: float | None
    radius_um: float
    loo_rows: list[dict[str, object]]


def prepare_analysis_table(data_root: Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map({"Female": 0.0, "Male": 1.0}).astype(float)

    rows: list[dict[str, object]] = []
    for row in cohort.itertuples(index=False):
        slide_name = str(row.slide_name)
        donor_id = str(row.donor_id)
        slide_path = data_root / slide_name
        svs_path = data_root / slide_name.replace(".svs.zarr", ".svs")
        features = extract_donor_features(slide_path=slide_path, svs_path=svs_path)
        rows.append({"donor_id": donor_id, "slide_name": slide_name, **features})

    features_df = pd.DataFrame(rows)
    return cohort.merge(features_df, on=["donor_id", "slide_name"], how="left")


def loo_predictions(df: pd.DataFrame, *, feature_col: str) -> tuple[float | None, list[dict[str, object]]]:
    needed = ["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]
    frame = df.loc[:, needed].dropna().reset_index(drop=True)
    if len(frame) < 3:
        return None, []

    predicted_residuals: list[float] = []
    actual_residuals: list[float] = []
    rows: list[dict[str, object]] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train_conf = _design_matrix(train[CONFOUND_COLUMNS])
        x_test_conf = _design_matrix(test[CONFOUND_COLUMNS])

        beta_feature, *_ = np.linalg.lstsq(x_train_conf, train[feature_col].to_numpy(dtype=float), rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train_conf, train[OUTCOME_COLUMN].to_numpy(dtype=float), rcond=None)

        resid_feature_train = train[feature_col].to_numpy(dtype=float) - x_train_conf @ beta_feature
        resid_outcome_train = train[OUTCOME_COLUMN].to_numpy(dtype=float) - x_train_conf @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            pred_resid = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            resid_feature_test = test[feature_col].to_numpy(dtype=float) - x_test_conf @ beta_feature
            pred_resid = float(slope * resid_feature_test[0])

        actual_resid = float(test[OUTCOME_COLUMN].to_numpy(dtype=float)[0] - (x_test_conf @ beta_outcome)[0])
        pred_outcome = float((x_test_conf @ beta_outcome)[0] + pred_resid) if math.isfinite(pred_resid) else float("nan")

        predicted_residuals.append(pred_resid)
        actual_residuals.append(actual_resid)
        rows.append(
            {
                "donor_id": str(test.loc[0, "donor_id"]),
                "outcome": float(test.loc[0, OUTCOME_COLUMN]),
                "predicted": pred_outcome,
                "predicted_residual": pred_resid,
                "actual_residual": actual_resid,
                feature_col: float(test.loc[0, feature_col]),
            }
        )

    pred_arr = np.asarray(predicted_residuals, dtype=float)
    act_arr = np.asarray(actual_residuals, dtype=float)
    mask = np.isfinite(pred_arr) & np.isfinite(act_arr)
    if mask.sum() < 3 or np.std(pred_arr[mask]) == 0 or np.std(act_arr[mask]) == 0:
        loo_r = None
    else:
        loo_r = float(np.corrcoef(pred_arr[mask], act_arr[mask])[0, 1])
    return loo_r, rows


def evaluate_variation(df: pd.DataFrame, *, variation_name: str, feature_col: str) -> VariationEvaluation:
    clean = df.loc[:, ["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]].dropna().reset_index(drop=True)
    n_total = int(len(df))
    n_analyzable = int(len(clean))

    if n_analyzable < 3:
        return VariationEvaluation(
            name=variation_name,
            feature_column=feature_col,
            partial_r=None,
            p_value=None,
            n_analyzable=n_analyzable,
            n_total=n_total,
            selection_score=None,
            loo_predictive_r=None,
            is_loo_gap=None,
            gap_penalty=None,
            adjusted_score=None,
            radius_um=LOCAL_RADIUS_UM,
            loo_rows=[],
        )

    metric = partial_correlation(clean, feature_col=feature_col, outcome_col=OUTCOME_COLUMN, confounds=CONFOUND_COLUMNS)
    partial_r = _safe_float(metric.get("partial_r"))
    p_value = _safe_float(metric.get("p_value"))
    selection = _selection_score(partial_r, n_analyzable, n_total)

    loo_r, loo_rows = loo_predictions(df, feature_col=feature_col)
    gap, gap_penalty, adjusted = _adjusted_metrics(partial_r, loo_r)

    return VariationEvaluation(
        name=variation_name,
        feature_column=feature_col,
        partial_r=partial_r,
        p_value=p_value,
        n_analyzable=n_analyzable,
        n_total=n_total,
        selection_score=selection,
        loo_predictive_r=loo_r,
        is_loo_gap=gap,
        gap_penalty=gap_penalty,
        adjusted_score=adjusted,
        radius_um=LOCAL_RADIUS_UM,
        loo_rows=loo_rows,
    )


def _variation_sort_key(item: VariationEvaluation) -> tuple[float, float]:
    sel = item.selection_score if item.selection_score is not None and math.isfinite(item.selection_score) else float("-inf")
    loo = item.loo_predictive_r if item.loo_predictive_r is not None and math.isfinite(item.loo_predictive_r) else float("-inf")
    return (sel, loo)


def _write_report(
    *,
    df: pd.DataFrame,
    best: VariationEvaluation,
    ranked: list[VariationEvaluation],
    report_path: Path,
) -> None:
    loo_frame = pd.DataFrame(best.loo_rows)
    top_errors = []
    if not loo_frame.empty:
        loo_frame["abs_error"] = (loo_frame["predicted"] - loo_frame["outcome"]).abs()
        error_df = loo_frame.merge(
            df[["donor_id", "braak_numeric", "cerad_ordinal", "cognitive_status", "sex"]].drop_duplicates(),
            on="donor_id",
            how="left",
        ).sort_values("abs_error", ascending=False)
        top_errors = error_df.head(3).to_dict(orient="records")

    error_lines = []
    if top_errors:
        for row in top_errors:
            error_lines.append(
                f"- {row['donor_id']}: outcome {row['outcome']:.3f}, prediction {row['predicted']:.3f}, "
                f"{best.feature_column}={row[best.feature_column]:.3f}, Braak {row.get('braak_numeric')}, "
                f"CERAD {row.get('cerad_ordinal')}, status {row.get('cognitive_status')}."
            )
    else:
        error_lines.append("- No analyzable LOO error rows were available.")

    ranking_text = "\n".join(
        f"- {item.name}: partial_r={item.partial_r:.4f} selection={item.selection_score:.4f} loo_r={item.loo_predictive_r:.4f}"
        if item.partial_r is not None and item.selection_score is not None and item.loo_predictive_r is not None
        else f"- {item.name}: insufficient analyzable donors"
        for item in ranked
    )

    worked_why = (
        f"The winner, {best.name}, tracks the donor-level fraction of CA1 reactive astrocytes sitting in "
        f"pyramidal-neuron-poor microenvironments within a 35 µm CA1 niche. That is biologically aligned with a "
        f"reactive glial response to local neuron depletion or scar-like remodeling, so it matches the existing CA1 injury theme."
    )
    failed_why = ""
    if len(ranked) > 1:
        runner_up = ranked[1]
        failed_why = (
            f"The nearby alternative {runner_up.name} was weaker, suggesting that its threshold either labeled too many "
            f"reactive astrocytes as 'desert' cells or blurred the contrast between clearly neuron-depleted niches and more mixed neighborhoods."
        )
    else:
        failed_why = "No nearby alternative was available to compare."

    rationale = (
        f"{best.name} won because the threshold of {VARIATIONS[best.name]['threshold']:.0%} keeps the biomarker focused on "
        f"reactive astrocytes occupying sharply neuron-depleted CA1 pockets, instead of merely reactive-rich but still neuron-mixed neighborhoods. "
        f"That makes it a coherent local refinement of the accepted panel's broader CA1 neuron-loss and reactive-astrocyte features."
    )

    interpretation = (
        "The signal most plausibly represents a CA1 reactive astrocyte niche in which local pyramidal neurons have dropped out, "
        "leaving reactive astrocytes embedded in neuron-poor tissue. Population: CA1 reactive astrocytes. "
        "Niche: their 35 µm CA1 neighborhood. Summary: fraction of reactive centers below a local pyramidal-neuron fraction threshold. "
        "Observable pattern: clusters of reactive astrocytes occupying locally neuron-sparse CA1 pockets."
    )

    next_text = (
        "Next local sweep: keep the same CA1 reactive-astrocyte-centered niche, but test whether weighting each center by its local reactive-cell density "
        "or requiring a minimum neighbor count sharpens the signal among the donors with the largest LOO residuals."
    )

    report = f"""## Summary
One sentence: tested the {HYPOTHESIS_FAMILY} family, {best.name} won, and its local winner selection score was {best.selection_score:.4f}.

## Metrics
Winning variation `{best.name}` used feature column `{best.feature_column}` with IS partial r {best.partial_r:.4f}, selection score {best.selection_score:.4f}, LOO predictive r {best.loo_predictive_r:.4f}, IS-LOO gap {best.is_loo_gap:.4f}, gap penalty {best.gap_penalty:.4f}, and adjusted score {best.adjusted_score:.4f}.

Other tested variations:
{ranking_text}

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - {worked_why}
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - {failed_why}
3. Error pattern: which donors are consistently wrong and what they share
{chr(10).join(error_lines)}

## Rationale
{rationale}
If relevant to the panel, this niche-focused reactive-astrocyte measure appears more specific than a global reactive burden feature and could plausibly add information beyond the current panel if it is not fully redundant with prior CA1 reactive-purity and hypertrophy members.

## Interpretation
{interpretation}

## Next
{next_text}
"""
    report_path.write_text(report, encoding="utf-8")


def main() -> int:
    data_root = Path("/data")
    scratch = Path("/scratch")
    scratch.mkdir(parents=True, exist_ok=True)

    table_path = scratch / "donor_feature_table.csv"
    results_path = scratch / "results.json"
    report_path = scratch / "report.md"
    canonical_path = Path(__file__).with_name("canonical_variation.json")

    analysis = prepare_analysis_table(data_root)

    variation_results: list[VariationEvaluation] = []
    for variation_name, spec in VARIATIONS.items():
        variation_results.append(
            evaluate_variation(
                analysis,
                variation_name=variation_name,
                feature_col=str(spec["feature_column"]),
            )
        )

    ranked = sorted(variation_results, key=_variation_sort_key, reverse=True)
    best = ranked[0]

    canonical_payload = {
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": best.name,
        "feature_column": best.feature_column,
        "local_radius_um": LOCAL_RADIUS_UM,
        "threshold": float(VARIATIONS[best.name]["threshold"]),
    }
    canonical_path.write_text(json.dumps(canonical_payload, indent=2) + "\n", encoding="utf-8")

    extra_columns = [
        item.feature_column for item in ranked
    ] + ["ca1_cell_count", "ca1_reactive_count", "ca1_pyramidal_count", "valid_center_count", "radius_px"]
    write_donor_feature_table(
        table_path,
        analysis,
        feature_column=best.feature_column,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUND_COLUMNS,
        id_columns=["donor_id", "slide_name"],
        extra_columns=extra_columns,
    )

    ranked_payload = []
    for item in ranked:
        ranked_payload.append(
            {
                "name": item.name,
                "feature_column": item.feature_column,
                "partial_r": item.partial_r,
                "selection_score": item.selection_score,
                "loo_predictive_r": item.loo_predictive_r,
                "is_loo_gap": item.is_loo_gap,
                "gap": item.is_loo_gap,
                "penalty": item.gap_penalty,
                "gap_penalty": item.gap_penalty,
                "adjusted_score": item.adjusted_score,
                "p_value": item.p_value,
                "n_analyzable": item.n_analyzable,
                "n_total": item.n_total,
                "radius_um": item.radius_um,
                "threshold": float(VARIATIONS[item.name]["threshold"]),
            }
        )

    payload = {
        "status": "ok",
        "feature_name": FEATURE_NAME,
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": best.name,
        "feature_column": best.feature_column,
        "outcome": OUTCOME_COLUMN,
        "n_total": best.n_total,
        "n_analyzable": best.n_analyzable,
        "partial_r": best.partial_r,
        "p_value": best.p_value,
        "selection_score": best.selection_score,
        "loo_predictive_r": best.loo_predictive_r,
        "is_loo_gap": best.is_loo_gap,
        "gap_penalty": best.gap_penalty,
        "penalty": best.gap_penalty,
        "adjusted_score": best.adjusted_score,
        "ranked_variations": ranked_payload,
        "donor_feature_table": str(table_path),
        "loo_rows": best.loo_rows,
        "artifacts": {
            "donor_feature_table": str(table_path),
            "canonical_variation": str(canonical_path),
        },
    }
    results_path.write_text(json.dumps(payload, indent=2, default=_json_default) + "\n", encoding="utf-8")

    _write_report(df=analysis, best=best, ranked=ranked, report_path=report_path)

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best.name}")
    print(f"  IS partial r:      {float(best.partial_r):.4f}" if best.partial_r is not None else "  IS partial r:      nan")
    print(f"  Selection score:   {float(best.selection_score):.4f}" if best.selection_score is not None else "  Selection score:   nan")
    print(f"  LOO predictive r:  {float(best.loo_predictive_r):.4f}  (diagnostic)" if best.loo_predictive_r is not None else "  LOO predictive r:  nan  (diagnostic)")
    if best.is_loo_gap is None or best.gap_penalty is None:
        print("  IS-LOO Gap:        nan  (penalty=nan)")
    else:
        print(f"  IS-LOO Gap:        {best.is_loo_gap:.4f}  (penalty={best.gap_penalty:.4f})")
    print(f"  Adjusted Score:    {float(best.adjusted_score):.4f}" if best.adjusted_score is not None else "  Adjusted Score:    nan")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked:
        pr = item.partial_r if item.partial_r is not None else float("nan")
        ss = item.selection_score if item.selection_score is not None else float("nan")
        lr = item.loo_predictive_r if item.loo_predictive_r is not None else float("nan")
        print(f"  {item.name}  {pr:.4f}  {ss:.4f}  {lr:.4f}")
    print()
    print("PER-DONOR (LOO):")
    for row in best.loo_rows:
        value = row.get(best.feature_column, float("nan"))
        print(f"  {row['donor_id']}  {row['outcome']:.4f}  {row['predicted']:.4f}  {value:.4f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
PY
python /scratch/result.py
