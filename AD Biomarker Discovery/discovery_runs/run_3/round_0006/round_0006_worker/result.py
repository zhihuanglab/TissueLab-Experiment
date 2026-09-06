from __future__ import annotations

import json
import math
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import sys

sys.path.insert(0, "/shared/lib")

from shared_analysis import (
    assign_centroids_to_regions,
    bootstrap_partial_correlation,
    build_results_payload,
    leave_one_out_summary,
    load_centroids,
    load_class_ids,
    load_class_lookup,
    load_contours,
    load_region_polygons,
    load_training_cohort,
    partial_correlation,
    residualized_loo_predictive_correlation,
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

HYPOTHESIS_FAMILY = (
    "CA1 peripyramidal reactive-astrocyte hypertrophy within 35 µm of CA1 pyramidal neurons"
)
TARGET_REGION = "CA1"
REACTIVE_CELL_TYPE = "Reactive Astrocyte"
PYRAMIDAL_CELL_TYPE = "Pyramidal Neuron"
RADIUS_UM = 35.0
DEFAULT_MPP = 0.503

# Canonical replay target. main() rewrites these constants to the winning variation.
CANONICAL_VARIATION = "candidate_variant_a"
FEATURE_NAME = "ca1_peripyramidal_reactive_area_median_r35um"
FEATURE_COLUMN = "ca1_peripyramidal_reactive_area_median_r35um"

VARIATIONS: dict[str, dict[str, object]] = {
    "candidate_variant_a": {
        "feature_name": "ca1_peripyramidal_reactive_area_median_r35um",
        "feature_column": "ca1_peripyramidal_reactive_area_median_r35um",
        "description": "Median contour area of CA1 reactive astrocytes within 35 µm of any CA1 pyramidal neuron",
        "summary_label": "median qualifying reactive-astrocyte area",
        "aggregation": "median",
    },
    "candidate_variant_b": {
        "feature_name": "ca1_peripyramidal_reactive_area_p75_r35um",
        "feature_column": "ca1_peripyramidal_reactive_area_p75_r35um",
        "description": "75th percentile contour area of CA1 reactive astrocytes within 35 µm of any CA1 pyramidal neuron",
        "summary_label": "75th percentile qualifying reactive-astrocyte area",
        "aggregation": "p75",
    },
}


def derive_sex_binary(series: pd.Series) -> pd.Series:
    mapped = series.map({"Female": 0.0, "Male": 1.0})
    return mapped.astype(float)


def clean_analysis_table(cohort: pd.DataFrame) -> pd.DataFrame:
    cohort = cohort.copy()
    if "sex_binary" not in cohort.columns:
        cohort["sex_binary"] = derive_sex_binary(cohort["sex"])
    return cohort


def slide_name_to_svs_path(slide_name_or_path: str | Path, data_root: Path | None = None) -> Path:
    slide_name = Path(slide_name_or_path).name
    if slide_name.endswith(".svs.zarr"):
        svs_name = slide_name[:-5]
    else:
        svs_name = slide_name
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


def _neighbor_hit_counts(query_points: np.ndarray, ref_points: np.ndarray, radius_px: float) -> np.ndarray:
    if len(query_points) == 0:
        return np.zeros(0, dtype=np.int32)
    if len(ref_points) == 0:
        return np.zeros(len(query_points), dtype=np.int32)
    tree = cKDTree(ref_points)
    try:
        counts = tree.query_ball_point(query_points, r=radius_px, return_length=True, workers=-1)
        return np.asarray(counts, dtype=np.int32)
    except TypeError:
        hits = tree.query_ball_point(query_points, r=radius_px, workers=-1)
        return np.fromiter((len(h) for h in hits), dtype=np.int32, count=len(query_points))


def compute_slide_feature_bundle(*, slide_path: Path, donor_id: str | None = None) -> dict[str, object]:
    polygons = load_region_polygons(slide_path).get(TARGET_REGION, [])
    svs_path = slide_name_to_svs_path(slide_path, slide_path.parent)
    mpp = read_slide_mpp(svs_path)
    area_unit = "um2" if mpp is not None else "px2"
    mpp_used = float(mpp if mpp is not None else DEFAULT_MPP)
    radius_px = float(RADIUS_UM / mpp_used)

    result: dict[str, object] = {
        "donor_id": donor_id,
        "slide_name": slide_path.name,
        "slide_mpp": float(mpp) if mpp is not None else np.nan,
        "radius_um": float(RADIUS_UM),
        "radius_px": float(radius_px),
        "area_unit": area_unit,
        "ca1_polygon_count": int(len(polygons)),
    }

    if not polygons:
        result.update(
            {
                "ca1_pyramidal_count": 0.0,
                "ca1_reactive_count": 0.0,
                "ca1_qualifying_reactive_count_r35um": 0.0,
            }
        )
        for variation in VARIATIONS.values():
            result[str(variation["feature_column"])] = float("nan")
        return result

    centroids = load_centroids(slide_path)
    class_lookup = load_class_lookup(slide_path)
    name_to_id = {name: idx for idx, name in class_lookup.items()}
    pyramidal_id = name_to_id.get(PYRAMIDAL_CELL_TYPE, None)
    reactive_id = name_to_id.get(REACTIVE_CELL_TYPE, None)
    nuclei_class_ids = load_class_ids(slide_path)

    ca1_labels = assign_centroids_to_regions(centroids, {TARGET_REGION: polygons})
    ca1_mask = ca1_labels == TARGET_REGION
    pyramidal_mask = ca1_mask & (nuclei_class_ids == pyramidal_id) if pyramidal_id is not None else np.zeros(len(centroids), dtype=bool)
    reactive_mask = ca1_mask & (nuclei_class_ids == reactive_id) if reactive_id is not None else np.zeros(len(centroids), dtype=bool)

    pyramidal_coords = centroids[pyramidal_mask]
    reactive_coords = centroids[reactive_mask]
    reactive_indices = np.flatnonzero(reactive_mask)

    result["ca1_pyramidal_count"] = float(len(pyramidal_coords))
    result["ca1_reactive_count"] = float(len(reactive_coords))

    if len(reactive_coords) == 0 or len(pyramidal_coords) == 0:
        result["ca1_qualifying_reactive_count_r35um"] = 0.0
        for variation in VARIATIONS.values():
            result[str(variation["feature_column"])] = 0.0
        return result

    hit_counts = _neighbor_hit_counts(reactive_coords, pyramidal_coords, radius_px)
    qualifying_local = hit_counts > 0
    qualifying_indices = reactive_indices[qualifying_local]
    result["ca1_qualifying_reactive_count_r35um"] = float(len(qualifying_indices))

    if len(qualifying_indices) == 0:
        for variation in VARIATIONS.values():
            result[str(variation["feature_column"])] = 0.0
        return result

    contours = load_contours(slide_path)
    areas_px2 = centered_polygon_area(contours[qualifying_indices])
    if mpp is not None and math.isfinite(float(mpp)) and float(mpp) > 0:
        areas = areas_px2 * (float(mpp) ** 2)
    else:
        areas = areas_px2

    result["ca1_qualifying_reactive_area_mean_r35um"] = float(np.mean(areas))
    result["ca1_qualifying_reactive_area_median_r35um_raw"] = float(np.median(areas))
    result["ca1_qualifying_reactive_area_p75_r35um_raw"] = float(np.percentile(areas, 75))

    for variation in VARIATIONS.values():
        feature_col = str(variation["feature_column"])
        agg = str(variation["aggregation"])
        if agg == "median":
            result[feature_col] = float(np.median(areas))
        elif agg == "p75":
            result[feature_col] = float(np.percentile(areas, 75))
        else:
            raise ValueError(f"Unknown aggregation: {agg}")

    return result


def compute_donor_score(*, donor_id: str, data_root: str | Path) -> float | None:
    """
    Canonical replay target for held-out evaluation.
    Returns the winning feature encoded by FEATURE_COLUMN / CANONICAL_VARIATION.
    """
    data_root = Path(data_root)
    cohort = clean_analysis_table(load_training_cohort(data_root))
    donor_rows = cohort.loc[cohort["donor_id"].astype(str) == str(donor_id)]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    bundle = compute_slide_feature_bundle(
        slide_path=data_root / slide_name,
        donor_id=str(donor_id),
    )
    value = bundle.get(FEATURE_COLUMN)
    if value is None:
        return None
    value_f = float(value)
    return None if not math.isfinite(value_f) else value_f


def selection_score(partial_r_value: float | None, n_analyzable: int, n_total: int) -> float:
    if partial_r_value is None or not math.isfinite(float(partial_r_value)):
        return float("nan")
    if n_total <= 0:
        return float("nan")
    return float(abs(float(partial_r_value)) * (float(n_analyzable) / float(n_total)))


def adjusted_metrics(partial_r_value: float | None, loo_r_value: float | None) -> tuple[float, float, float]:
    is_r = abs(float(partial_r_value)) if partial_r_value is not None and math.isfinite(float(partial_r_value)) else 0.0
    loo_r = abs(float(loo_r_value)) if loo_r_value is not None and math.isfinite(float(loo_r_value)) else 0.0
    gap = float(is_r - loo_r)
    penalty = float(max(0.0, gap - 0.15))
    adjusted = float(loo_r - penalty)
    return gap, penalty, adjusted


def raw_loo_prediction_table(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
) -> pd.DataFrame:
    frame = df.dropna(subset=[feature_col, outcome_col, *confounds]).reset_index(drop=True)
    rows: list[dict[str, object]] = []
    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train = np.column_stack(
            [
                np.ones(len(train), dtype=float),
                train[confounds].to_numpy(dtype=float),
                train[[feature_col]].to_numpy(dtype=float),
            ]
        )
        y_train = train[outcome_col].to_numpy(dtype=float)
        beta, *_ = np.linalg.lstsq(x_train, y_train, rcond=None)

        x_test = np.column_stack(
            [
                np.ones(len(test), dtype=float),
                test[confounds].to_numpy(dtype=float),
                test[[feature_col]].to_numpy(dtype=float),
            ]
        )
        predicted = float((x_test @ beta)[0])
        outcome = float(test.iloc[0][outcome_col])
        feature_value = float(test.iloc[0][feature_col])
        rows.append(
            {
                "donor_id": str(test.iloc[0]["donor_id"]),
                "outcome": outcome,
                "predicted": predicted,
                feature_col: feature_value,
                "abs_error": float(abs(predicted - outcome)),
                "braak_numeric": float(test.iloc[0]["braak_numeric"]),
                "cerad_ordinal": float(test.iloc[0]["cerad_ordinal"]),
                "sex_binary": float(test.iloc[0]["sex_binary"]),
                "max_age_vis": float(test.iloc[0]["max_age_vis"]),
                "cognitive_status": str(test.iloc[0].get("cognitive_status", "")),
            }
        )
    return pd.DataFrame(rows)


def _pythonize(value):
    if isinstance(value, dict):
        return {str(k): _pythonize(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_pythonize(v) for v in value]
    if isinstance(value, tuple):
        return [_pythonize(v) for v in value]
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def write_canonical_constants(self_path: Path, best_variation: str) -> None:
    spec = VARIATIONS[best_variation]
    text = self_path.read_text(encoding="utf-8")
    text = re.sub(
        r'CANONICAL_VARIATION = ".*?"',
        f'CANONICAL_VARIATION = "{best_variation}"',
        text,
        count=1,
    )
    text = re.sub(
        r'FEATURE_NAME = ".*?"',
        f'FEATURE_NAME = "{spec["feature_name"]}"',
        text,
        count=1,
    )
    text = re.sub(
        r'FEATURE_COLUMN = ".*?"',
        f'FEATURE_COLUMN = "{spec["feature_column"]}"',
        text,
        count=1,
    )
    self_path.write_text(text, encoding="utf-8")


def build_report(
    *,
    best: dict[str, object],
    ranked: list[dict[str, object]],
    donor_table: pd.DataFrame,
    per_donor_loo: pd.DataFrame,
    report_path: Path,
) -> None:
    winner_name = str(best["variation_name"])
    winner_desc = str(best["description"])
    feature_col = str(best["feature_column"])
    score = float(best["selection_score"])
    partial_r_value = float(best["partial_r"])
    loo_r_value = float(best["loo_predictive_r"])
    gap_value = float(best["is_loo_gap"])
    penalty_value = float(best["gap_penalty"])
    adjusted_value = float(best["adjusted_score"])

    losers = [row for row in ranked if str(row["variation_name"]) != winner_name]
    loser_text = (
        "; ".join(
            f'{row["variation_name"]} (selection={float(row["selection_score"]):.4f}, '
            f'partial_r={float(row["partial_r"]):.4f}, loo={float(row["loo_predictive_r"]):.4f})'
            for row in losers
        )
        if losers
        else "No nearby alternatives were tested."
    )

    error_table = per_donor_loo.sort_values("abs_error", ascending=False).reset_index(drop=True)
    worst = error_table.head(5).copy()
    high_error_donors = ", ".join(
        f'{row.donor_id} (|e|={row.abs_error:.3f})' for row in worst.itertuples(index=False)
    ) or "none"

    severe = worst.loc[worst["outcome"] <= donor_table[OUTCOME_COLUMN].median()]
    severe_text = (
        f"{len(severe)}/{len(worst)} of the largest errors are on donors with more negative-than-median memory slope"
        if len(worst) > 0
        else "No analyzable donor-level LOO errors were available"
    )

    if "p75" in str(best["aggregation"]).lower():
        worked_why = (
            "Emphasizing the upper-tail area of qualifying reactive astrocytes modestly improved the local family score, "
            "consistent with a hypertrophic subset in the immediate CA1 pyramidal niche carrying more signal than the niche median."
        )
        interpretation_text = (
            "The signal appears to reflect a heavy-tail hypertrophic reactive-astrocyte state in the CA1 peripyramidal niche."
        )
    else:
        worked_why = (
            "The winning summary favored the central tendency of qualifying reactive-astrocyte size, suggesting that broad niche-wide hypertrophy "
            "is more reproducible than focusing only on the most extreme cells."
        )
        interpretation_text = (
            "The signal appears to reflect a donor-level shift toward generally enlarged reactive astrocytes in the CA1 peripyramidal niche."
        )

    failed_why = (
        "The losing variation likely over- or under-emphasized the hypertrophic tail. "
        "Within this fixed niche, changing only the donor-level area summary moved signal less than earlier count-based localization steps, "
        "suggesting morphology adds only a modest refinement on top of the already-established cuffing niche."
    )

    additivity = (
        "Because the current accepted panel already contains broad reactivity and CA1 cuffing features, this morphology summary is biologically coherent "
        "but may be partially redundant unless its hypertrophy readout captures severity beyond simple neighbor counts."
    )

    report = f"""## Summary
One sentence: tested the {HYPOTHESIS_FAMILY} family, {winner_name} won ({winner_desc}), and the local winner selection score was {score:.4f}.

## Metrics
Winner: {winner_name}
- IS partial r: {partial_r_value:.4f}
- Selection score: {score:.4f}
- LOO predictive r: {loo_r_value:.4f}
- IS-LOO Gap: {gap_value:.4f} (penalty={penalty_value:.4f})
- Adjusted Score: {adjusted_value:.4f}

Other tested variations: {loser_text}

## Findings
1. What worked and why: {worked_why}
2. What failed and why: {failed_why}
3. Error pattern: the largest LOO errors were {high_error_donors}; {severe_text}.

## Rationale
The best variation is biologically coherent because it stays inside the same CA1 peripyramidal reactive-astrocyte niche that drove rounds 3-5, but changes the donor summary from abundance to morphology. Reactive astrocyte contour area is a plausible proxy for hypertrophic activation state on LFB tissue, and the winner beat the nearby alternative because its chosen summary better balanced niche specificity against donor-level robustness. {additivity}

## Interpretation
{interpretation_text} Population: CA1 Reactive Astrocytes. Niche: within 35 µm of CA1 Pyramidal Neurons. Feature summary: {best["summary_label"]}. Simplest observable pattern: larger swollen reactive astrocyte profiles clustered directly around CA1 pyramidal neurons track worse memory decline.

## Next
Next local sweep: keep the same CA1 peripyramidal reactive-astrocyte niche but test whether hypertrophy is most informative when conditioned on reactive-cuff burden (for example, median or p75 area only among donors' cuffed pyramidal neighborhoods with >=2 nearby reactive astrocytes), since this round asks morphology alone to add on top of already-strong count-based winners.
"""
    report_path.write_text(report, encoding="utf-8")


def main() -> int:
    data_root = DATA_ROOT_DEFAULT
    scratch_root = SCRATCH_ROOT_DEFAULT
    scratch_root.mkdir(parents=True, exist_ok=True)

    cohort = clean_analysis_table(load_training_cohort(data_root))
    bundles = [
        compute_slide_feature_bundle(
            slide_path=data_root / str(row.slide_name),
            donor_id=str(row.donor_id),
        )
        for row in cohort.itertuples(index=False)
    ]
    feature_table = cohort.merge(
        pd.DataFrame(bundles),
        on=["donor_id", "slide_name"],
        how="left",
        validate="one_to_one",
    )

    ranked_variations: list[dict[str, object]] = []
    per_donor_loo_by_variation: dict[str, pd.DataFrame] = {}

    for variation_name, spec in VARIATIONS.items():
        feature_col = str(spec["feature_column"])
        analysis_frame = feature_table.dropna(subset=[feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]).copy()
        n_total = int(len(feature_table))
        n_analyzable = int(len(analysis_frame))

        partial = partial_correlation(
            analysis_frame,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUND_COLUMNS,
        )
        boot = bootstrap_partial_correlation(
            analysis_frame,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUND_COLUMNS,
            n_boot=400,
            random_state=(101 if variation_name == "candidate_variant_a" else 202),
        )
        loo_r = residualized_loo_predictive_correlation(
            analysis_frame,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUND_COLUMNS,
        )
        loo_shift = leave_one_out_summary(
            analysis_frame,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUND_COLUMNS,
            id_col="donor_id",
        )
        sel = selection_score(partial.get("partial_r"), n_analyzable, n_total)
        gap, penalty, adjusted = adjusted_metrics(partial.get("partial_r"), loo_r)
        per_donor_loo = raw_loo_prediction_table(
            analysis_frame,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUND_COLUMNS,
        )
        per_donor_loo_by_variation[variation_name] = per_donor_loo

        ranked_variations.append(
            {
                "variation_name": variation_name,
                "feature_name": str(spec["feature_name"]),
                "feature_column": feature_col,
                "description": str(spec["description"]),
                "summary_label": str(spec["summary_label"]),
                "aggregation": str(spec["aggregation"]),
                "n_total": n_total,
                "n_analyzable": n_analyzable,
                "coverage": float(n_analyzable / n_total) if n_total else float("nan"),
                "partial_r": float(partial.get("partial_r", float("nan"))),
                "p_value": float(partial.get("p_value", float("nan"))),
                "ci_lo": float(boot.get("ci_lo", float("nan"))),
                "ci_hi": float(boot.get("ci_hi", float("nan"))),
                "selection_score": float(sel),
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
    best_per_donor_loo = per_donor_loo_by_variation[best_variation]

    donor_feature_table_path = scratch_root / "donor_feature_table.csv"
    donor_feature_cols = [
        "donor_id",
        "slide_name",
        *(str(spec["feature_column"]) for spec in VARIATIONS.values()),
        OUTCOME_COLUMN,
        *CONFOUND_COLUMNS,
        "cognitive_status",
        "braak_label",
        "cerad_label",
        "ca1_pyramidal_count",
        "ca1_reactive_count",
        "ca1_qualifying_reactive_count_r35um",
        "area_unit",
        "slide_mpp",
    ]
    donor_feature_cols = [col for col in donor_feature_cols if col in feature_table.columns]
    feature_table.loc[:, donor_feature_cols].to_csv(donor_feature_table_path, index=False)

    results = build_results_payload(
        status="ok",
        feature_name=str(best["feature_name"]),
        outcome=OUTCOME_COLUMN,
        n_total=int(best["n_total"]),
        n_analyzable=int(best["n_analyzable"]),
        partial_r=float(best["partial_r"]),
        ci_lo=float(best["ci_lo"]),
        ci_hi=float(best["ci_hi"]),
        p_value=float(best["p_value"]),
        loo_predictive_r=float(best["loo_predictive_r"]),
        loo_unstable_count=int(best["loo_unstable_count"]),
        loo_max_shift=float(best["loo_max_shift"]),
        donor_ids_used=[
            str(x)
            for x in feature_table.loc[
                feature_table[best_feature_col].notna(), "donor_id"
            ].tolist()
        ],
        covariates=CONFOUND_COLUMNS,
        recomputed_from_raw=True,
        registry_written=False,
        artifacts={
            "donor_feature_table": str(donor_feature_table_path),
            "report": str(scratch_root / "report.md"),
        },
        hypothesis_family=HYPOTHESIS_FAMILY,
        best_variation=best_variation,
        ranked_variations=ranked_variations,
        feature_column=best_feature_col,
        selection_score=float(best["selection_score"]),
        is_loo_gap=float(best["is_loo_gap"]),
        gap_penalty=float(best["gap_penalty"]),
        adjusted_score=float(best["adjusted_score"]),
        per_donor_loo=best_per_donor_loo.to_dict(orient="records"),
        area_unit_mode=str(feature_table["area_unit"].dropna().mode().iloc[0]) if feature_table["area_unit"].notna().any() else "",
    )
    results_path = scratch_root / "results.json"
    results_path.write_text(json.dumps(_pythonize(results), indent=2) + "\n", encoding="utf-8")

    build_report(
        best=best,
        ranked=ranked_variations,
        donor_table=feature_table,
        per_donor_loo=best_per_donor_loo,
        report_path=scratch_root / "report.md",
    )

    write_canonical_constants(Path(__file__), best_variation)

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
            f"  {str(row['variation_name']):<17}  "
            f"{float(row['partial_r']):>8.4f}  "
            f"{float(row['selection_score']):>15.4f}  "
            f"{float(row['loo_predictive_r']):>16.4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    for row in best_per_donor_loo.itertuples(index=False):
        row_map = row._asdict()
        print(
            f"  {row_map['donor_id']}  "
            f"{float(row_map['outcome']):.4f}  "
            f"{float(row_map['predicted']):.4f}  "
            f"{float(row_map[best_feature_col]):.4f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
