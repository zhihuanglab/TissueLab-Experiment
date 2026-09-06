from __future__ import annotations

import json
import math
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import zarr
from matplotlib.path import Path as MplPath
from scipy.stats import pearsonr

warnings.filterwarnings("ignore")

HYPOTHESIS_FAMILY = "CA1-localized pyramidal neuron fraction"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
MIN_CLASSIFIED_CELLS = 100

KNOWN_REGIONS = ("CA1", "CA2", "CA3", "CA4", "DG", "EC", "SB", "TEC")
UNKNOWN_LABELS = {"", "unknown", "unlabeled", "unlabelled"}

VARIATIONS = [
    {
        "name": "candidate_variant_a",
        "description": "CA1 pyramidal neuron fraction",
        "regions": ("CA1",),
        "feature_column": "ca1_pyramidal_fraction",
    },
    {
        "name": "candidate_variant_b",
        "description": "CA1_CA2 pyramidal neuron fraction",
        "regions": ("CA1", "CA2"),
        "feature_column": "ca1_ca2_pyramidal_fraction",
    },
]

CANONICAL_VARIATION = "candidate_variant_a"
FEATURE_NAME = "ca1_pyramidal_neuron_fraction"
FEATURE_COLUMN = "ca1_pyramidal_fraction"


def _is_ignored_name(name: str) -> bool:
    return name.startswith("._") or name.startswith(".DS_Store")


def _decode_scalar(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if hasattr(value, "item"):
        try:
            return _decode_scalar(value.item())
        except Exception:
            return value
    return value


def slide_id_from_name(name: str) -> str:
    base = Path(name).name
    for suffix in (".svs.zarr", ".svs"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def canonicalize_region_label(label: str | None) -> str:
    token = (label or "").strip().upper()
    for region in sorted(KNOWN_REGIONS, key=len, reverse=True):
        if token.startswith(region):
            return region
    return token or "UNKNOWN"


def load_training_cohort(data_root: str | Path) -> pd.DataFrame:
    data_root = Path(data_root)
    cohort_path = data_root / "training_cohort.csv"
    cohort = pd.read_csv(cohort_path)
    if "slide_id" not in cohort.columns and "slide_name" in cohort.columns:
        cohort["slide_id"] = cohort["slide_name"].map(slide_id_from_name)
    return cohort


def open_slide_zarr(zarr_path: str | Path):
    return zarr.open(str(Path(zarr_path)), mode="r")


def load_centroids(zarr_path: str | Path) -> np.ndarray:
    root = open_slide_zarr(zarr_path)
    return np.asarray(root["SegmentationNode"]["centroids"][:], dtype=np.float32)


def load_class_ids(zarr_path: str | Path) -> np.ndarray:
    root = open_slide_zarr(zarr_path)
    return np.asarray(root["ClassificationNode"]["nuclei_class_id"][:], dtype=np.int32)


def load_class_names(zarr_path: str | Path) -> list[str]:
    root = open_slide_zarr(zarr_path)
    return [str(_decode_scalar(x)) for x in root["ClassificationNode"]["nuclei_class_name"][:]]


def _extract_polygon_points(annotation: dict[str, Any], *, scale: float) -> np.ndarray:
    geometry = annotation.get("target", {}).get("selector", {}).get("geometry", {})
    points = geometry.get("points")
    if points is None and geometry.get("coordinates"):
        coordinates = geometry["coordinates"]
        if coordinates and coordinates[0]:
            points = coordinates[0]
    if not points:
        return np.empty((0, 2), dtype=np.float32)
    polygon = np.asarray(points, dtype=np.float32)
    if polygon.ndim != 2 or polygon.shape[1] != 2:
        return np.empty((0, 2), dtype=np.float32)
    return polygon / float(scale)


def load_region_polygons(zarr_path: str | Path, *, scale: float = 16.0) -> dict[str, list[np.ndarray]]:
    root = open_slide_zarr(zarr_path)
    annotations = root["CustomAnnotations"]
    grouped: dict[str, list[np.ndarray]] = {}
    for key in sorted(annotations.group_keys()):
        if _is_ignored_name(key):
            continue
        group = annotations[key]
        raw_json = _decode_scalar(group["annotation_json"][()])
        if not raw_json:
            continue
        try:
            annotation = json.loads(raw_json)
        except json.JSONDecodeError:
            continue
        raw_label = str(_decode_scalar(group["comment"][()]) or "").strip()
        polygon = _extract_polygon_points(annotation, scale=scale)
        if len(polygon) < 3:
            continue
        grouped.setdefault(canonicalize_region_label(raw_label), []).append(polygon)
    return grouped


def assign_centroids_to_regions(
    centroids: np.ndarray,
    region_polygons: dict[str, list[np.ndarray]],
) -> np.ndarray:
    centroids = np.asarray(centroids, dtype=np.float32)
    labels = np.full(len(centroids), "", dtype=object)
    for region, polygons in region_polygons.items():
        region_mask = np.zeros(len(centroids), dtype=bool)
        for polygon in polygons:
            if len(polygon) < 3:
                continue
            region_mask |= MplPath(polygon, closed=True).contains_points(centroids)
        assignable = region_mask & (labels == "")
        labels[assignable] = region
    labels[labels == ""] = None
    return labels


def _class_id_sets(zarr_path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    class_names = load_class_names(zarr_path)
    valid_ids: list[int] = []
    pyramidal_ids: list[int] = []
    for idx, name in enumerate(class_names):
        token = str(name).strip().lower()
        if token not in UNKNOWN_LABELS:
            valid_ids.append(idx)
        if token == "pyramidal neuron":
            pyramidal_ids.append(idx)
    return np.asarray(valid_ids, dtype=np.int32), np.asarray(pyramidal_ids, dtype=np.int32)


def extract_donor_region_counts(zarr_path: str | Path) -> dict[str, int]:
    centroids = load_centroids(zarr_path)
    class_ids = load_class_ids(zarr_path)
    valid_ids, pyramidal_ids = _class_id_sets(zarr_path)
    region_labels = assign_centroids_to_regions(centroids, load_region_polygons(zarr_path, scale=16.0))

    valid_mask = np.isin(class_ids, valid_ids)
    pyramidal_mask = np.isin(class_ids, pyramidal_ids)

    counts: dict[str, int] = {}
    for region in ("CA1", "CA2"):
        region_mask = region_labels == region
        denom = int(np.sum(region_mask & valid_mask))
        numer = int(np.sum(region_mask & valid_mask & pyramidal_mask))
        counts[f"{region.lower()}_classified"] = denom
        counts[f"{region.lower()}_pyramidal"] = numer
        counts[f"{region.lower()}_all_cells"] = int(np.sum(region_mask))
    return counts


def compute_fraction(counts: dict[str, int], regions: tuple[str, ...], *, min_cells: int = MIN_CLASSIFIED_CELLS) -> float:
    denom = int(sum(counts.get(f"{r.lower()}_classified", 0) for r in regions))
    numer = int(sum(counts.get(f"{r.lower()}_pyramidal", 0) for r in regions))
    if denom < min_cells:
        return float("nan")
    return float(numer / denom)


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
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    counts = extract_donor_region_counts(slide_path)
    regions = next(v["regions"] for v in VARIATIONS if v["name"] == CANONICAL_VARIATION)
    score = compute_fraction(counts, regions)
    return None if not np.isfinite(score) else float(score)


def design_matrix(confounds: pd.DataFrame) -> np.ndarray:
    mat = confounds.astype(float).to_numpy()
    intercept = np.ones((len(confounds), 1), dtype=float)
    return np.hstack([intercept, mat])


def residualize(values: pd.Series | np.ndarray, confounds: pd.DataFrame) -> np.ndarray:
    y = np.asarray(values, dtype=float)
    x = design_matrix(confounds)
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ beta


def partial_correlation(df: pd.DataFrame, *, feature_col: str, outcome_col: str, confounds: list[str]) -> dict[str, float]:
    frame = df[[feature_col, outcome_col, *confounds]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(frame) < 3:
        return {"n": float(len(frame)), "partial_r": float("nan"), "p_value": float("nan")}
    feat_resid = residualize(frame[feature_col], frame[confounds])
    out_resid = residualize(frame[outcome_col], frame[confounds])
    if np.std(feat_resid) == 0 or np.std(out_resid) == 0:
        return {"n": float(len(frame)), "partial_r": float("nan"), "p_value": float("nan")}
    corr, p_value = pearsonr(feat_resid, out_resid)
    return {"n": float(len(frame)), "partial_r": float(corr), "p_value": float(p_value)}


def leave_one_out_predictions(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
) -> tuple[float, pd.DataFrame]:
    frame = df[["donor_id", feature_col, outcome_col, *confounds]].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    rows: list[dict[str, float | str]] = []
    pred_resids: list[float] = []
    actual_resids: list[float] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train = design_matrix(train[confounds])
        x_test = design_matrix(test[confounds])

        beta_feature, *_ = np.linalg.lstsq(x_train, train[feature_col].to_numpy(dtype=float), rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train, train[outcome_col].to_numpy(dtype=float), rcond=None)

        resid_feature_train = train[feature_col].to_numpy(dtype=float) - x_train @ beta_feature
        resid_outcome_train = train[outcome_col].to_numpy(dtype=float) - x_train @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        resid_feature_test = test[feature_col].to_numpy(dtype=float) - x_test @ beta_feature
        actual_resid_test = test[outcome_col].to_numpy(dtype=float) - x_test @ beta_outcome
        if denom <= 0:
            pred_resid_test = np.array([0.0], dtype=float)
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            pred_resid_test = slope * resid_feature_test

        raw_pred = float((x_test @ beta_outcome)[0] + pred_resid_test[0])

        pred_resids.append(float(pred_resid_test[0]))
        actual_resids.append(float(actual_resid_test[0]))
        rows.append(
            {
                "donor_id": str(test.loc[0, "donor_id"]),
                "outcome": float(test.loc[0, outcome_col]),
                "predicted": raw_pred,
                "actual_residual": float(actual_resid_test[0]),
                "predicted_residual": float(pred_resid_test[0]),
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

    per_donor = pd.DataFrame(rows)
    per_donor["error"] = per_donor["predicted"] - per_donor["outcome"]
    per_donor["abs_error"] = per_donor["error"].abs()
    return loo_r, per_donor


def leave_one_out_influence(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
) -> dict[str, Any]:
    frame = df[["donor_id", feature_col, outcome_col, *confounds]].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    base = partial_correlation(frame, feature_col=feature_col, outcome_col=outcome_col, confounds=confounds)
    if frame.empty or not np.isfinite(base["partial_r"]):
        return {"max_shift": float("nan"), "unstable_count": 0, "unstable_donors": []}
    rows = []
    for idx in range(len(frame)):
        reduced = frame.drop(index=idx).reset_index(drop=True)
        result = partial_correlation(reduced, feature_col=feature_col, outcome_col=outcome_col, confounds=confounds)
        delta = float(result["partial_r"] - base["partial_r"]) if np.isfinite(result["partial_r"]) else float("nan")
        rows.append({"donor_id": frame.loc[idx, "donor_id"], "partial_r": result["partial_r"], "delta_from_full": delta})
    loo = pd.DataFrame(rows)
    shifts = loo["delta_from_full"].abs()
    unstable = loo.loc[shifts > 0.10].copy()
    return {
        "max_shift": float(shifts.max()) if len(shifts) else float("nan"),
        "unstable_count": int((shifts > 0.10).sum()) if len(shifts) else 0,
        "unstable_donors": unstable.to_dict(orient="records"),
    }


def evaluate_variation(
    df: pd.DataFrame,
    variation: dict[str, Any],
    *,
    outcome_col: str = OUTCOME_COLUMN,
    confounds: list[str] = CONFOUND_COLUMNS,
) -> dict[str, Any]:
    feature_col = variation["feature_column"]
    stats = partial_correlation(df, feature_col=feature_col, outcome_col=outcome_col, confounds=confounds)
    n_analyzable = int(stats["n"])
    n_total = int(len(df))
    coverage = float(n_analyzable / n_total) if n_total else float("nan")
    selection_score = abs(float(stats["partial_r"])) * coverage if np.isfinite(stats["partial_r"]) else float("nan")
    loo_r, per_donor = leave_one_out_predictions(df, feature_col=feature_col, outcome_col=outcome_col, confounds=confounds)
    is_loo_gap = abs(float(stats["partial_r"])) - abs(float(loo_r)) if np.isfinite(stats["partial_r"]) and np.isfinite(loo_r) else float("nan")
    penalty = max(0.0, is_loo_gap) if np.isfinite(is_loo_gap) else float("nan")
    adjusted_score = selection_score - penalty if np.isfinite(selection_score) and np.isfinite(penalty) else float("nan")
    influence = leave_one_out_influence(df, feature_col=feature_col, outcome_col=outcome_col, confounds=confounds)
    return {
        "variation_name": variation["name"],
        "description": variation["description"],
        "feature_column": feature_col,
        "regions": list(variation["regions"]),
        "n_analyzable": n_analyzable,
        "n_total": n_total,
        "coverage": coverage,
        "partial_r": float(stats["partial_r"]),
        "p_value": float(stats["p_value"]),
        "selection_score": float(selection_score),
        "loo_predictive_r": float(loo_r),
        "is_loo_gap": float(is_loo_gap),
        "penalty": float(penalty),
        "adjusted_score": float(adjusted_score),
        "max_shift": influence["max_shift"],
        "unstable_count": influence["unstable_count"],
        "unstable_donors": influence["unstable_donors"],
        "per_donor_loo": per_donor.to_dict(orient="records"),
    }


def build_feature_table(data_root: str | Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map(lambda x: 1.0 if str(x).strip().lower().startswith("m") else 0.0)
    feature_rows: list[dict[str, Any]] = []
    for _, row in cohort.iterrows():
        slide_path = Path(data_root) / str(row["slide_name"])
        counts = extract_donor_region_counts(slide_path)
        feature_rows.append(
            {
                "donor_id": row["donor_id"],
                "slide_name": row["slide_name"],
                **counts,
                "ca1_pyramidal_fraction": compute_fraction(counts, ("CA1",)),
                "ca1_ca2_pyramidal_fraction": compute_fraction(counts, ("CA1", "CA2")),
            }
        )
    features = pd.DataFrame(feature_rows)
    merged = cohort.merge(features, on=["donor_id", "slide_name"], how="left")
    return merged


def emit_stdout(best: dict[str, Any], ranked: list[dict[str, Any]]) -> None:
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
    for row in ranked:
        print(
            f"  {row['variation_name']}  "
            f"{row['partial_r']:.4f}  "
            f"{row['selection_score']:.4f}  "
            f"{row['loo_predictive_r']:.4f}"
        )
    print()
    feature_col = best["feature_column"]
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {feature_col}")
    for row in best["per_donor_loo"]:
        print(
            f"  {row['donor_id']}  "
            f"{row['outcome']:.4f}  "
            f"{row['predicted']:.4f}  "
            f"{row[feature_col]:.6f}"
        )


def main() -> None:
    data_root = Path("/data")
    out_root = Path("/scratch")

    donor_table = build_feature_table(data_root)
    donor_table.to_csv(out_root / "donor_feature_table.csv", index=False)

    ranked = [evaluate_variation(donor_table, variation) for variation in VARIATIONS]
    ranked.sort(
        key=lambda row: (
            -np.inf if not np.isfinite(row["selection_score"]) else -row["selection_score"],
            -np.inf if not np.isfinite(row["partial_r"]) else -abs(row["partial_r"]),
        )
    )
    best = ranked[0]

    results_payload = {
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": best["variation_name"],
        "best_variation_description": best["description"],
        "feature_column": best["feature_column"],
        "partial_r": best["partial_r"],
        "selection_score": best["selection_score"],
        "loo_predictive_r": best["loo_predictive_r"],
        "is_loo_gap": best["is_loo_gap"],
        "penalty": best["penalty"],
        "adjusted_score": best["adjusted_score"],
        "ranked_variations": ranked,
    }
    with open(out_root / "results.json", "w") as f:
        json.dump(results_payload, f, indent=2)

    emit_stdout(best, ranked)


if __name__ == "__main__":
    main()
