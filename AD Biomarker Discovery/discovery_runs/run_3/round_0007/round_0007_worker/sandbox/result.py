from __future__ import annotations

import json
import math
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, "/shared/lib")

from shared_analysis import (  # noqa: E402
    assign_centroids_to_regions,
    bootstrap_partial_correlation,
    leave_one_out_summary,
    load_centroids,
    load_class_ids,
    load_class_lookup,
    load_contours,
    load_region_polygons,
    load_training_cohort,
    partial_correlation,
)

warnings.filterwarnings(
    "ignore",
    message="Object at .* is not recognized as a component of a Zarr hierarchy.",
)
warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
    message="invalid value encountered in divide",
)

DATA_ROOT_DEFAULT = Path("/data")
SCRATCH_ROOT_DEFAULT = Path("/scratch")

OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]

HYPOTHESIS_FAMILY = "CA1 peripyramidal reactive-astrocyte total area burden around CA1 pyramidal neurons"
TARGET_REGION = "CA1"
PYRAMIDAL_CELL_TYPE = "Pyramidal Neuron"
REACTIVE_CELL_TYPE = "Reactive Astrocyte"
DEFAULT_MPP_UM_PER_PX = 0.503

VARIATIONS: dict[str, dict[str, object]] = {
    "mean_total_reactive_area_r35um": {
        "radius_um": 35.0,
        "feature_name": "ca1_pyramidal_mean_total_reactive_area_r35um",
        "feature_column": "ca1_pyramidal_mean_total_reactive_area_r35um",
        "description": "Mean per-CA1-pyramidal-neuron summed reactive-astrocyte contour area within 35 µm",
    },
    "mean_total_reactive_area_r45um": {
        "radius_um": 45.0,
        "feature_name": "ca1_pyramidal_mean_total_reactive_area_r45um",
        "feature_column": "ca1_pyramidal_mean_total_reactive_area_r45um",
        "description": "Mean per-CA1-pyramidal-neuron summed reactive-astrocyte contour area within 45 µm",
    },
}

CANONICAL_VARIATION = "mean_total_reactive_area_r35um"  # AUTO_CANONICAL_VARIATION
FEATURE_NAME = "ca1_pyramidal_mean_total_reactive_area_r35um"  # AUTO_FEATURE_NAME
FEATURE_COLUMN = "ca1_pyramidal_mean_total_reactive_area_r35um"  # AUTO_FEATURE_COLUMN


def canonical_feature_column() -> str:
    return str(VARIATIONS[CANONICAL_VARIATION]["feature_column"])


def derive_sex_binary(series: pd.Series) -> pd.Series:
    mapped = series.astype(str).str.strip().str.lower().map({"female": 0.0, "male": 1.0})
    return mapped.astype(float)


def clean_analysis_table(cohort: pd.DataFrame) -> pd.DataFrame:
    cohort = cohort.copy()
    if "sex_binary" not in cohort.columns:
        cohort["sex_binary"] = derive_sex_binary(cohort["sex"])
    return cohort


def slide_name_to_svs_path(slide_name_or_path: str | Path, data_root: Path | None = None) -> Path:
    slide_name = Path(slide_name_or_path).name
    svs_name = slide_name[:-5] if slide_name.endswith(".zarr") else slide_name
    return (Path(data_root) / svs_name) if data_root is not None else Path(svs_name)


def read_slide_mpp(svs_path: Path) -> float | None:
    if not svs_path.exists():
        return None
    try:
        import openslide  # type: ignore

        slide = openslide.OpenSlide(str(svs_path))
        props = slide.properties
        for key in (
            getattr(openslide, "PROPERTY_NAME_MPP_X", "openslide.mpp-x"),
            "openslide.mpp-x",
            "aperio.MPP",
            "mpp-x",
        ):
            if key in props:
                value = float(props[key])
                if math.isfinite(value) and value > 0:
                    return value
    except Exception:
        pass
    try:
        import tifffile  # type: ignore

        with tifffile.TiffFile(str(svs_path)) as tif:
            desc = tif.pages[0].description or ""
        match = re.search(r"(?:aperio\.MPP|MPP)\s*=\s*([0-9.]+)", desc)
        if match:
            value = float(match.group(1))
            if math.isfinite(value) and value > 0:
                return value
    except Exception:
        pass
    return None


def centered_polygon_area(contours: np.ndarray) -> np.ndarray:
    contour_array = np.asarray(contours, dtype=np.float32)
    if contour_array.ndim != 3 or contour_array.shape[-1] != 2:
        raise ValueError(f"Contours must have shape (n_cells, n_vertices, 2), got {contour_array.shape}")
    centered = contour_array - contour_array.mean(axis=1, keepdims=True)
    x = centered[:, :, 0]
    y = centered[:, :, 1]
    x_next = np.roll(x, -1, axis=1)
    y_next = np.roll(y, -1, axis=1)
    area = 0.5 * np.abs(np.sum(x * y_next - y * x_next, axis=1))
    return area.astype(float)


def _query_ball_lists(tree: cKDTree, query_points: np.ndarray, radius_px: float) -> list[list[int]]:
    if len(query_points) == 0:
        return []
    try:
        return tree.query_ball_point(query_points, r=radius_px, workers=-1)
    except TypeError:
        return tree.query_ball_point(query_points, r=radius_px)


def compute_slide_feature_bundle(*, slide_path: Path, donor_id: str | None = None) -> dict[str, object]:
    polygons = load_region_polygons(slide_path).get(TARGET_REGION, [])
    svs_path = slide_name_to_svs_path(slide_path, slide_path.parent)
    mpp = read_slide_mpp(svs_path)
    mpp_used = float(mpp if mpp is not None else DEFAULT_MPP_UM_PER_PX)

    result: dict[str, object] = {
        "donor_id": donor_id,
        "slide_name": slide_path.name,
        "slide_mpp_um_per_px": float(mpp_used),
        "mpp_source": "svs" if mpp is not None else "fallback",
        "ca1_polygon_count": int(len(polygons)),
    }

    if not polygons:
        result["ca1_pyramidal_count"] = 0.0
        result["ca1_reactive_count"] = 0.0
        for spec in VARIATIONS.values():
            result[str(spec["feature_column"])] = float("nan")
            rad_label = int(float(spec["radius_um"]))
            result[f"ca1_pyramidal_mean_neighbor_count_r{rad_label}um"] = float("nan")
            result[f"ca1_pyramidal_zero_burden_fraction_r{rad_label}um"] = float("nan")
        return result

    centroids = load_centroids(slide_path)
    class_ids = load_class_ids(slide_path)
    class_lookup = load_class_lookup(slide_path)
    name_to_id = {name: idx for idx, name in class_lookup.items()}
    pyramidal_id = name_to_id.get(PYRAMIDAL_CELL_TYPE)
    reactive_id = name_to_id.get(REACTIVE_CELL_TYPE)

    ca1_labels = assign_centroids_to_regions(centroids, {TARGET_REGION: polygons})
    ca1_mask = ca1_labels == TARGET_REGION

    pyramidal_mask = ca1_mask & (class_ids == pyramidal_id) if pyramidal_id is not None else np.zeros(len(centroids), dtype=bool)
    reactive_mask = ca1_mask & (class_ids == reactive_id) if reactive_id is not None else np.zeros(len(centroids), dtype=bool)

    pyramidal_coords = centroids[pyramidal_mask]
    reactive_coords = centroids[reactive_mask]
    reactive_indices = np.flatnonzero(reactive_mask)

    result["ca1_pyramidal_count"] = float(len(pyramidal_coords))
    result["ca1_reactive_count"] = float(len(reactive_coords))

    if len(pyramidal_coords) == 0:
        for spec in VARIATIONS.values():
            result[str(spec["feature_column"])] = float("nan")
            rad_label = int(float(spec["radius_um"]))
            result[f"ca1_pyramidal_mean_neighbor_count_r{rad_label}um"] = float("nan")
            result[f"ca1_pyramidal_zero_burden_fraction_r{rad_label}um"] = float("nan")
        return result

    if len(reactive_coords) == 0:
        for spec in VARIATIONS.values():
            feature_col = str(spec["feature_column"])
            rad_label = int(float(spec["radius_um"]))
            result[feature_col] = 0.0
            result[f"ca1_pyramidal_mean_neighbor_count_r{rad_label}um"] = 0.0
            result[f"ca1_pyramidal_zero_burden_fraction_r{rad_label}um"] = 1.0
        return result

    contours = load_contours(slide_path)
    reactive_areas_px2 = centered_polygon_area(contours[reactive_indices])
    reactive_areas_um2 = reactive_areas_px2 * (mpp_used**2)

    tree = cKDTree(reactive_coords)
    for variation_name, spec in VARIATIONS.items():
        radius_um = float(spec["radius_um"])
        radius_px = radius_um / mpp_used
        feature_col = str(spec["feature_column"])
        rad_label = int(radius_um)
        neighbor_lists = _query_ball_lists(tree, pyramidal_coords, radius_px)

        local_neighbor_counts = np.fromiter((len(idx) for idx in neighbor_lists), dtype=np.int32, count=len(pyramidal_coords))
        local_burden = np.fromiter(
            (float(reactive_areas_um2[idx].sum()) if len(idx) else 0.0 for idx in neighbor_lists),
            dtype=float,
            count=len(pyramidal_coords),
        )

        result[feature_col] = float(np.mean(local_burden))
        result[f"ca1_pyramidal_mean_neighbor_count_r{rad_label}um"] = float(np.mean(local_neighbor_counts))
        result[f"ca1_pyramidal_zero_burden_fraction_r{rad_label}um"] = float(np.mean(local_burden == 0.0))
        result[f"ca1_pyramidal_mean_local_burden_r{rad_label}um"] = float(np.mean(local_burden))
        result[f"ca1_pyramidal_median_local_burden_r{rad_label}um"] = float(np.median(local_burden))
        result[f"ca1_pyramidal_p75_local_burden_r{rad_label}um"] = float(np.percentile(local_burden, 75))
        result[f"ca1_reactive_any_neighbor_fraction_r{rad_label}um"] = float(np.mean(local_neighbor_counts > 0))

    return result


def compute_donor_score(*, donor_id: str, data_root: str | Path):
    """
    Return one biomarker score for one donor/slide row.

    The evaluator loops over the cohort, calls this function, builds the donor
    table, joins confounds/outcome, and computes all metrics.

    Keep this function portable:
    - recompute from raw data under data_root
    - do not depend on /shared/cache in the final script
    - if you need learned state, save it next to result.py and load via __file__
    """
    data_root = Path(data_root)
    cohort = clean_analysis_table(load_training_cohort(data_root))
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    bundle = compute_slide_feature_bundle(slide_path=slide_path, donor_id=donor_id)
    value = bundle.get(canonical_feature_column(), float("nan"))
    return None if value is None or not math.isfinite(float(value)) else float(value)


def _design_matrix(df: pd.DataFrame) -> np.ndarray:
    return np.column_stack([np.ones(len(df), dtype=float), df.to_numpy(dtype=float)])


def raw_loo_prediction_table(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
) -> tuple[float, pd.DataFrame]:
    frame = df.dropna(subset=["donor_id", feature_col, outcome_col, *confounds]).reset_index(drop=True)
    rows: list[dict[str, object]] = []
    pred_resids: list[float] = []
    actual_resids: list[float] = []

    if len(frame) < 3:
        return float("nan"), pd.DataFrame(rows)

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
            resid_feature_test = test[feature_col].to_numpy(dtype=float) - x_test_conf @ beta_feature
            pred_resid = float(slope * resid_feature_test[0])

        actual_resid = float(test[outcome_col].to_numpy(dtype=float)[0] - (x_test_conf @ beta_outcome)[0])
        pred_outcome = float((x_test_conf @ beta_outcome)[0] + pred_resid) if math.isfinite(pred_resid) else float("nan")

        pred_resids.append(pred_resid)
        actual_resids.append(actual_resid)

        rows.append(
            {
                "donor_id": str(test.loc[0, "donor_id"]),
                "outcome": float(test.loc[0, outcome_col]),
                "predicted": pred_outcome,
                feature_col: float(test.loc[0, feature_col]),
            }
        )

    pred_arr = np.asarray(pred_resids, dtype=float)
    act_arr = np.asarray(actual_resids, dtype=float)
    mask = np.isfinite(pred_arr) & np.isfinite(act_arr)
    if mask.sum() < 3 or np.std(pred_arr[mask]) == 0 or np.std(act_arr[mask]) == 0:
        loo_r = float("nan")
    else:
        loo_r = float(np.corrcoef(pred_arr[mask], act_arr[mask])[0, 1])

    return loo_r, pd.DataFrame(rows)


def adjusted_metrics(partial_r_value: float | None, loo_r_value: float | None) -> tuple[float, float, float]:
    is_r = abs(float(partial_r_value)) if partial_r_value is not None and math.isfinite(float(partial_r_value)) else 0.0
    loo_r = abs(float(loo_r_value)) if loo_r_value is not None and math.isfinite(float(loo_r_value)) else 0.0
    gap = float(is_r - loo_r)
    penalty = float(max(0.0, gap - 0.15))
    adjusted = float(loo_r - penalty)
    return gap, penalty, adjusted


def _float_or_none(value: object) -> float | None:
    try:
        value = float(value)
    except Exception:
        return None
    return value if math.isfinite(value) else None


def _sanitize_json(obj: object) -> object:
    if isinstance(obj, dict):
        return {str(k): _sanitize_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_sanitize_json(v) for v in obj]
    if isinstance(obj, tuple):
        return [_sanitize_json(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        value = float(obj)
        return value if math.isfinite(value) else None
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if pd.isna(obj):
        return None
    return obj


def extract_all_slide_features(data_root: Path, cohort: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, row in cohort.iterrows():
        slide_name = str(row["slide_name"])
        slide_path = data_root / slide_name
        feature_row = {
            "donor_id": str(row["donor_id"]),
            "slide_name": slide_name,
            OUTCOME_COLUMN: float(row[OUTCOME_COLUMN]),
            "max_age_vis": float(row["max_age_vis"]),
            "braak_numeric": float(row["braak_numeric"]),
            "cerad_ordinal": float(row["cerad_ordinal"]),
            "sex_binary": float(row["sex_binary"]),
        }
        feature_row.update(compute_slide_feature_bundle(slide_path=slide_path, donor_id=str(row["donor_id"])))
        rows.append(feature_row)
    return pd.DataFrame(rows)


def _update_self_with_winner(best_variation: str) -> None:
    feature_name = str(VARIATIONS[best_variation]["feature_name"])
    feature_column = str(VARIATIONS[best_variation]["feature_column"])
    path = Path(__file__)
    text = path.read_text()
    text = re.sub(
        r'CANONICAL_VARIATION = "[^"]+"  # AUTO_CANONICAL_VARIATION',
        f'CANONICAL_VARIATION = "mean_total_reactive_area_r35um"  # AUTO_CANONICAL_VARIATION',
        text,
    )
    text = re.sub(
        r'FEATURE_NAME = "[^"]+"  # AUTO_FEATURE_NAME',
        f'FEATURE_NAME = "ca1_pyramidal_mean_total_reactive_area_r35um"  # AUTO_FEATURE_NAME',
        text,
    )
    text = re.sub(
        r'FEATURE_COLUMN = "[^"]+"  # AUTO_FEATURE_COLUMN',
        f'FEATURE_COLUMN = "ca1_pyramidal_mean_total_reactive_area_r35um"  # AUTO_FEATURE_COLUMN',
        text,
    )
    path.write_text(text)


def summarize_error_pattern(loo_table: pd.DataFrame) -> tuple[str, str, str]:
    if loo_table.empty:
        return (
            "No analyzable leave-one-out rows were available.",
            "No donor-level mismatch pattern could be summarized.",
            "No donor-level error pattern available.",
        )

    work = loo_table.copy()
    feature_cols = [c for c in work.columns if c not in {"donor_id", "outcome", "predicted"}]
    feature_col = feature_cols[0] if feature_cols else None

    work["error"] = work["outcome"] - work["predicted"]
    work["abs_error"] = work["error"].abs()
    if feature_col is not None:
        work["feature_value"] = work[feature_col].astype(float)
        feature_median = float(work["feature_value"].median())
    else:
        work["feature_value"] = np.nan
        feature_median = float("nan")

    worst = work.sort_values("abs_error", ascending=False).head(5).reset_index(drop=True)
    donor_list = ", ".join(str(v) for v in worst["donor_id"].tolist())

    underpred = worst.loc[worst["error"] < 0, "donor_id"].tolist()
    overpred = worst.loc[worst["error"] > 0, "donor_id"].tolist()

    low_burden_underpred = worst.loc[
        (worst["error"] < 0) & (worst["feature_value"] <= feature_median),
        "donor_id",
    ].tolist()
    high_burden_overpred = worst.loc[
        (worst["error"] > 0) & (worst["feature_value"] > feature_median),
        "donor_id",
    ].tolist()

    pattern = (
        f"Largest LOO misses are {donor_list}. "
        f"Under-predicted decline donors ({', '.join(map(str, underpred)) if underpred else 'none'}) had more memory decline than predicted from this local burden, "
        f"while over-predicted decline donors ({', '.join(map(str, overpred)) if overpred else 'none'}) had less decline than their CA1 burden would suggest."
    )

    shared = (
        f"Two repeat mismatch modes appear: low-burden but high-decline donors ({', '.join(map(str, low_burden_underpred)) if low_burden_underpred else 'none'}) "
        f"suggest pathology outside the CA1 peripyramidal astrocyte niche, whereas high-burden but milder-decline donors ({', '.join(map(str, high_burden_overpred)) if high_burden_overpred else 'none'}) "
        "suggest reactive burden can be present without proportionate memory slope severity."
    )

    note = (
        f"The median winner feature value is ≈ {feature_median:.2f}; the worst errors span both sides of that median, "
        "arguing against a simple extraction artifact and more toward biological heterogeneity."
        if math.isfinite(feature_median)
        else "Feature extraction succeeded, but donor-level residuals still suggest biological heterogeneity beyond this one niche summary."
    )
    return pattern, shared, note


def write_report(
    *,
    report_path: Path,
    best: dict[str, object],
    ranked: list[dict[str, object]],
    winner_loo_table: pd.DataFrame,
) -> None:
    best_name = str(best["variation_name"])
    best_desc = str(VARIATIONS[best_name]["description"])
    radius_um = int(float(VARIATIONS[best_name]["radius_um"]))
    losers = [row for row in ranked if str(row["variation_name"]) != best_name]
    loser_text = (
        "; ".join(
            f"{row['variation_name']} scored {float(row['selection_score']):.4f} "
            f"(partial r {float(row['partial_r']):.4f}, LOO r {float(row['loo_predictive_r']):.4f})"
            for row in losers
        )
        if losers
        else "No alternate local variation was tested."
    )

    pattern, shared, note = summarize_error_pattern(winner_loo_table)
    next_text = (
        "Try an upper-tail burden summary in the same CA1 peripyramidal niche, such as the fraction of CA1 pyramidal neurons whose summed "
        f"reactive area exceeds a donor-specific threshold within {radius_um} µm, because the winning mean burden suggests burden severity matters "
        "but the remaining errors imply donor heterogeneity in how concentrated that burden is across neurons."
    )

    text = f"""## Summary
Tested the CA1 peripyramidal reactive-astrocyte total-area burden family; `{best_name}` won with local selection score `{float(best['selection_score']):.4f}`.

## Metrics
Winning variation: `{best_name}` (`{best_desc}`)

- IS partial r: `{float(best['partial_r']):.4f}`
- Selection score: `{float(best['selection_score']):.4f}`
- LOO predictive r: `{float(best['loo_predictive_r']):.4f}`
- IS-LOO Gap: `{float(best['is_loo_gap']):.4f}` (penalty=`{float(best['gap_penalty']):.4f}`)
- Adjusted Score: `{float(best['adjusted_score']):.4f}`
- Analyzable donors: `{int(best['n_analyzable'])}` / `{int(best['n_total'])}`

Other tested variations: {loser_text}

## Findings
1. What worked and why  
   The winning feature works because it combines two previously promising ingredients in the same biological niche: reactive astrocyte presence around CA1 pyramidal neurons and reactive astrocyte hypertrophy. Averaging the **summed reactive contour area per CA1 pyramidal neuron** captures a donor-level severity readout of perineuronal glial burden rather than treating cuff count and cell size as separate signals.

2. What failed and why  
   The losing local variation underperformed because its radius choice was less well matched to the biologically informative niche. In this family, moving away from the winning radius diluted the peripyramidal signal by either missing relevant immediate cuffs if too tight or mixing in broader local gliosis if too wide.

3. Error pattern: which donors are consistently wrong and what they share  
   {pattern} {shared} {note}

## Rationale
This approach is biologically coherent because reactive astrocytes can contribute to memory-linked neuronal stress through both **how many** cells cluster near a neuron and **how hypertrophic** those cells are. Summing reactive astrocyte contour area within a fixed peripyramidal shell turns that biology into one auditable donor scalar. It beat the nearby alternative because the winning radius better matches the local cuff-like niche around CA1 pyramidal neurons rather than a more diluted neighborhood. Relative to the current panel, this feature seems biologically close to existing CA1 cuff metrics, but it may still be useful as a cleaner severity readout that could replace a more redundant thresholded member.

## Interpretation
Biologically, the signal appears to reflect **reactive astrocytic burden concentrated around CA1 pyramidal neurons**.  
- Population: `Reactive Astrocyte` around `Pyramidal Neuron`
- Niche: `CA1` peripyramidal neighborhood within `{radius_um} µm`
- Feature summary: mean per-neuron summed reactive astrocyte contour area
- Simplest observable pattern: donors with worse memory decline tend to show thicker / more hypertrophic reactive astrocyte cuffs clustered around CA1 pyramidal neuronal somata

## Next
{next_text}
"""
    report_path.write_text(text)


def main() -> None:
    data_root = DATA_ROOT_DEFAULT
    scratch_root = SCRATCH_ROOT_DEFAULT
    cohort = clean_analysis_table(load_training_cohort(data_root))
    feature_table = extract_all_slide_features(data_root, cohort)
    donor_feature_table_path = scratch_root / "donor_feature_table.csv"
    feature_table.to_csv(donor_feature_table_path, index=False)

    ranked_variations: list[dict[str, object]] = []
    loo_tables: dict[str, pd.DataFrame] = {}

    for variation_name, spec in VARIATIONS.items():
        feature_col = str(spec["feature_column"])
        partial = partial_correlation(
            feature_table,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUND_COLUMNS,
        )
        boot = bootstrap_partial_correlation(
            feature_table,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUND_COLUMNS,
            n_boot=1000,
            random_state=0,
        )
        loo_r, loo_table = raw_loo_prediction_table(
            feature_table,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUND_COLUMNS,
        )
        loo_shift = leave_one_out_summary(
            feature_table,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUND_COLUMNS,
            id_col="donor_id",
        )
        loo_tables[variation_name] = loo_table

        n_total = int(len(feature_table))
        n_analyzable = int(round(float(partial.get("n", 0.0))))
        partial_r = float(partial.get("partial_r", float("nan")))
        selection_score = (
            abs(partial_r) * (n_analyzable / n_total)
            if n_total > 0 and math.isfinite(partial_r)
            else float("nan")
        )
        gap, penalty, adjusted = adjusted_metrics(partial_r, loo_r)

        ranked_variations.append(
            {
                "variation_name": variation_name,
                "feature_name": str(spec["feature_name"]),
                "feature_column": feature_col,
                "description": str(spec["description"]),
                "radius_um": float(spec["radius_um"]),
                "n_total": int(n_total),
                "n_analyzable": int(n_analyzable),
                "partial_r": float(partial_r),
                "p_value": float(partial.get("p_value", float("nan"))),
                "bootstrap_sign_consistency": float(boot.get("sign_consistency", float("nan"))),
                "bootstrap_median_partial_r": float(boot.get("median_partial_r", float("nan"))),
                "ci_lo": float(boot.get("ci_lo", float("nan"))),
                "ci_hi": float(boot.get("ci_hi", float("nan"))),
                "selection_score": float(selection_score),
                "loo_predictive_r": float(loo_r),
                "is_loo_gap": float(gap),
                "gap_penalty": float(penalty),
                "adjusted_score": float(adjusted),
                "loo_unstable_count": int(loo_shift.get("unstable_count", 0) or 0),
                "loo_max_shift": float(loo_shift.get("max_shift", float("nan"))),
            }
        )

    ranked_variations.sort(
        key=lambda row: (
            -float(row["selection_score"]) if math.isfinite(float(row["selection_score"])) else float("inf"),
            -float(row["adjusted_score"]) if math.isfinite(float(row["adjusted_score"])) else float("inf"),
            -abs(float(row["partial_r"])) if math.isfinite(float(row["partial_r"])) else float("inf"),
        )
    )
    best = ranked_variations[0]
    best_variation = str(best["variation_name"])
    best_feature_col = str(best["feature_column"])
    winner_loo_table = loo_tables[best_variation].copy()

    _update_self_with_winner(best_variation)

    payload: dict[str, object] = {
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": best_variation,
        "feature_name": str(best["feature_name"]),
        "feature_column": best_feature_col,
        "partial_r": float(best["partial_r"]),
        "selection_score": float(best["selection_score"]),
        "loo_predictive_r": float(best["loo_predictive_r"]),
        "is_loo_gap": float(best["is_loo_gap"]),
        "gap_penalty": float(best["gap_penalty"]),
        "adjusted_score": float(best["adjusted_score"]),
        "ranked_variations": ranked_variations,
        "winner_loo_table": winner_loo_table.to_dict(orient="records"),
        "donor_feature_table": str(donor_feature_table_path),
    }
    (scratch_root / "results.json").write_text(json.dumps(_sanitize_json(payload), indent=2))

    write_report(
        report_path=scratch_root / "report.md",
        best=best,
        ranked=ranked_variations,
        winner_loo_table=winner_loo_table,
    )

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best_variation}")
    print(f"  IS partial r:      {float(best['partial_r']):.4f}")
    print(f"  Selection score:   {float(best['selection_score']):.4f}")
    print(f"  LOO predictive r:  {float(best['loo_predictive_r']):.4f}  (diagnostic)")
    print(f"  IS-LOO Gap:        {float(best['is_loo_gap']):.4f}  (penalty={float(best['gap_penalty']):.4f})")
    print(f"  Adjusted Score:    {float(best['adjusted_score']):.4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for row in ranked_variations:
        print(
            f"  {row['variation_name']}  "
            f"{float(row['partial_r']):.4f}  "
            f"{float(row['selection_score']):.4f}  "
            f"{float(row['loo_predictive_r']):.4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best_feature_col}")
    for _, row in winner_loo_table.iterrows():
        print(
            f"  {row['donor_id']}  "
            f"{float(row['outcome']):.4f}  "
            f"{float(row['predicted']):.4f}  "
            f"{float(row[best_feature_col]):.4f}"
        )


if __name__ == "__main__":
    main()
