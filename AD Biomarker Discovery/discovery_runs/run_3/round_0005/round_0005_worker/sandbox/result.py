from __future__ import annotations

import json
import math
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore", category=UserWarning)

HYPOTHESIS_FAMILY = "ca1-reactive-dominant-cuff"
CANONICAL_VARIATION = "candidate_variant_b"
FEATURE_NAME = "ca1_pyramidal_supermajority_reactive_cuff_fraction_r35um_ge2"
FEATURE_COLUMN = "feature__ca1_pyramidal_supermajority_reactive_cuff_fraction_r35um_ge2"

TARGET_REGION = "CA1"
CENTER_TYPE = "Pyramidal Neuron"
ASTRO_TYPE = "Astrocyte"
REACTIVE_ASTRO_TYPE = "Reactive Astrocyte"
RADIUS_UM = 35.0
MIN_TOTAL_ASTRO = 2

VARIATION_SPECS = {
    "candidate_variant_a": {
        "feature_name": "ca1_pyramidal_majority_reactive_cuff_fraction_r35um_ge2",
        "feature_column": "feature__ca1_pyramidal_majority_reactive_cuff_fraction_r35um_ge2",
        "description": "CA1 pyramidal neurons with >=2 nearby CA1 astrocytes within 35 um and reactive fraction > 0.50",
        "threshold": 0.50,
        "op": ">",
    },
    "candidate_variant_b": {
        "feature_name": "ca1_pyramidal_supermajority_reactive_cuff_fraction_r35um_ge2",
        "feature_column": "feature__ca1_pyramidal_supermajority_reactive_cuff_fraction_r35um_ge2",
        "description": "CA1 pyramidal neurons with >=2 nearby CA1 astrocytes within 35 um and reactive fraction >= 0.67",
        "threshold": 0.67,
        "op": ">=",
    },
}


def decode_class_names(name_array) -> list[str]:
    out = []
    for x in name_array[:]:
        if isinstance(x, (bytes, np.bytes_)):
            out.append(x.decode("utf-8"))
        else:
            out.append(str(x))
    return out


def load_slide_mpp(svs_path: Path) -> float:
    try:
        import openslide

        slide = openslide.OpenSlide(str(svs_path))
        px = slide.properties.get(openslide.PROPERTY_NAME_MPP_X)
        py = slide.properties.get(openslide.PROPERTY_NAME_MPP_Y)
        vals = [float(v) for v in [px, py] if v is not None]
        if vals:
            return float(np.mean(vals))
    except Exception:
        pass

    try:
        import tifffile

        with tifffile.TiffFile(str(svs_path)) as tif:
            desc = tif.pages[0].description or ""
        m = re.search(r"\bMPP\s*=\s*([0-9.]+)", desc)
        if m:
            return float(m.group(1))
    except Exception:
        pass

    return 0.5


def load_region_polygons(root: zarr.Group, region_name: str) -> list[np.ndarray]:
    polygons: list[np.ndarray] = []
    ann_root = root.get("CustomAnnotations")
    if ann_root is None:
        return polygons
    for key in ann_root.keys():
        if not str(key).startswith(f"{region_name}_"):
            continue
        try:
            annotation_json = ann_root[key]["annotation_json"][()]
            if not annotation_json:
                continue
            obj = json.loads(str(annotation_json))
            pts = obj["target"]["selector"]["geometry"]["points"]
            arr = np.asarray(pts, dtype=np.float64)
            if arr.ndim == 2 and arr.shape[1] == 2 and len(arr) >= 3:
                polygons.append(arr / 16.0)
        except Exception:
            continue
    return polygons


def points_in_any_polygon(points: np.ndarray, polygons: list[np.ndarray]) -> np.ndarray:
    if len(points) == 0 or not polygons:
        return np.zeros(len(points), dtype=bool)
    try:
        from matplotlib.path import Path as MplPath

        inside = np.zeros(len(points), dtype=bool)
        for poly in polygons:
            inside |= MplPath(poly).contains_points(points)
        return inside
    except Exception:
        from skimage.measure import points_in_poly

        inside = np.zeros(len(points), dtype=bool)
        for poly in polygons:
            inside |= points_in_poly(points, poly)
        return inside


def neighbor_count_per_point(query_coords: np.ndarray, ref_coords: np.ndarray, radius_px: float) -> np.ndarray:
    if len(query_coords) == 0:
        return np.zeros(0, dtype=np.int32)
    if len(ref_coords) == 0:
        return np.zeros(len(query_coords), dtype=np.int32)
    tree = cKDTree(ref_coords)
    try:
        counts = tree.query_ball_point(query_coords, r=radius_px, return_length=True)
        return np.asarray(counts, dtype=np.int32)
    except TypeError:
        hits = tree.query_ball_point(query_coords, r=radius_px)
        return np.fromiter((len(h) for h in hits), dtype=np.int32, count=len(query_coords))


def evaluate_variation(total_astro: np.ndarray, reactive_astro: np.ndarray, denom_n: int, *, threshold: float, op: str) -> float:
    if denom_n <= 0:
        return float("nan")
    eligible = total_astro >= MIN_TOTAL_ASTRO
    reactive_frac = np.zeros(denom_n, dtype=np.float64)
    nz = total_astro > 0
    reactive_frac[nz] = reactive_astro[nz] / total_astro[nz]
    if op == ">":
        meet = eligible & (reactive_frac > threshold)
    elif op == ">=":
        meet = eligible & (reactive_frac >= threshold)
    else:
        raise ValueError(op)
    return float(meet.mean())


def compute_slide_variations(slide_path: Path) -> dict:
    root = zarr.open_group(str(slide_path), mode="r")
    polygons = load_region_polygons(root, TARGET_REGION)

    centroids = np.asarray(root["SegmentationNode"]["centroids"][:], dtype=np.float32)
    if centroids.shape[1] != 2:
        raise ValueError(f"Unexpected centroid shape: {centroids.shape}")

    class_ids = np.asarray(root["ClassificationNode"]["nuclei_class_id"][:], dtype=np.int32)
    class_names = decode_class_names(root["ClassificationNode"]["nuclei_class_name"])
    name_to_id = {name: i for i, name in enumerate(class_names)}

    ca1_mask = points_in_any_polygon(centroids, polygons)
    ca1_centroids = centroids[ca1_mask]
    ca1_class_ids = class_ids[ca1_mask]

    pyr_id = name_to_id.get(CENTER_TYPE, -999)
    astro_id = name_to_id.get(ASTRO_TYPE, -999)
    reactive_id = name_to_id.get(REACTIVE_ASTRO_TYPE, -999)

    pyr_coords = ca1_centroids[ca1_class_ids == pyr_id]
    astro_coords = ca1_centroids[np.isin(ca1_class_ids, [astro_id, reactive_id])]
    reactive_coords = ca1_centroids[ca1_class_ids == reactive_id]

    n_pyr = int(len(pyr_coords))
    mpp = load_slide_mpp(slide_path.with_suffix(""))
    radius_px = float(RADIUS_UM / mpp)

    total_astro = neighbor_count_per_point(pyr_coords, astro_coords, radius_px)
    reactive_astro = neighbor_count_per_point(pyr_coords, reactive_coords, radius_px)

    result = {
        "n_ca1_cells": int(len(ca1_centroids)),
        "n_ca1_pyramidal": n_pyr,
        "n_ca1_astro_total": int(len(astro_coords)),
        "n_ca1_reactive_astro": int(len(reactive_coords)),
        "radius_um": RADIUS_UM,
        "mpp": float(mpp),
        "radius_px": float(radius_px),
        "scores": {},
    }
    for name, spec in VARIATION_SPECS.items():
        result["scores"][name] = evaluate_variation(
            total_astro,
            reactive_astro,
            n_pyr,
            threshold=spec["threshold"],
            op=spec["op"],
        )
    return result


def load_training_cohort(data_root: str | Path) -> pd.DataFrame:
    cohort = pd.read_csv(Path(data_root) / "training_cohort.csv")
    cohort["sex_binary"] = cohort["sex"].map({"Female": 0.0, "Male": 1.0})
    return cohort


def compute_donor_variations(*, donor_id: str, data_root: str | Path) -> dict[str, float]:
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return {k: float("nan") for k in VARIATION_SPECS}
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    slide_result = compute_slide_variations(slide_path)
    return slide_result["scores"]


def compute_donor_score(*, donor_id: str, data_root: str | Path):
    scores = compute_donor_variations(donor_id=donor_id, data_root=data_root)
    value = scores.get(CANONICAL_VARIATION, float("nan"))
    if value is None or not np.isfinite(value):
        return None
    return float(value)


def _safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    mask = np.isfinite(a) & np.isfinite(b)
    a = a[mask]
    b = b[mask]
    if len(a) < 3:
        return float("nan")
    a = a - a.mean()
    b = b - b.mean()
    sa = float(np.sqrt(np.sum(a * a)))
    sb = float(np.sqrt(np.sum(b * b)))
    if sa == 0.0 or sb == 0.0:
        return float("nan")
    return float(np.sum(a * b) / (sa * sb))


def _fit_ols(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta


def _residualize(y: np.ndarray, conf: np.ndarray) -> np.ndarray:
    X = np.column_stack([np.ones(len(y)), conf])
    beta = _fit_ols(X, y)
    return y - X @ beta


def partial_correlation(feature: np.ndarray, outcome: np.ndarray, conf: np.ndarray) -> float:
    rx = _residualize(feature, conf)
    ry = _residualize(outcome, conf)
    return _safe_corr(rx, ry)


def loo_diagnostic(feature: np.ndarray, outcome: np.ndarray, conf: np.ndarray):
    n = len(feature)
    pred_full = np.full(n, np.nan, dtype=np.float64)
    pred_resid = np.full(n, np.nan, dtype=np.float64)
    actual_resid = np.full(n, np.nan, dtype=np.float64)

    for i in range(n):
        tr = np.ones(n, dtype=bool)
        tr[i] = False

        X_conf_tr = np.column_stack([np.ones(tr.sum()), conf[tr]])
        X_conf_te = np.concatenate([[1.0], conf[i]])

        beta_conf = _fit_ols(X_conf_tr, outcome[tr])
        conf_pred_i = float(X_conf_te @ beta_conf)

        X_full_tr = np.column_stack([np.ones(tr.sum()), conf[tr], feature[tr]])
        X_full_te = np.concatenate([[1.0], conf[i], [feature[i]]])

        beta_full = _fit_ols(X_full_tr, outcome[tr])
        full_pred_i = float(X_full_te @ beta_full)

        pred_full[i] = full_pred_i
        pred_resid[i] = full_pred_i - conf_pred_i
        actual_resid[i] = outcome[i] - conf_pred_i

    loo_r = _safe_corr(pred_resid, actual_resid)
    return loo_r, pred_full, pred_resid, actual_resid


def analyze_variation(df: pd.DataFrame, variation_name: str) -> dict:
    spec = VARIATION_SPECS[variation_name]
    feature_col = spec["feature_column"]
    conf_cols = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
    req_cols = conf_cols + ["slope_zmem0", feature_col]
    mask = np.ones(len(df), dtype=bool)
    for col in req_cols:
        mask &= pd.to_numeric(df[col], errors="coerce").notna().to_numpy()

    dfa = df.loc[mask].copy()
    n_total = int(len(df))
    n_analyzable = int(len(dfa))
    if n_analyzable < 5:
        return {
            "variation_name": variation_name,
            "feature_name": spec["feature_name"],
            "feature_column": feature_col,
            "partial_r": float("nan"),
            "selection_score": float("nan"),
            "loo_predictive_r": float("nan"),
            "is_loo_gap": float("nan"),
            "penalty": float("nan"),
            "adjusted_score": float("nan"),
            "n_total": n_total,
            "n_analyzable": n_analyzable,
            "per_donor_loo": [],
        }

    x = dfa[feature_col].to_numpy(dtype=np.float64)
    y = dfa["slope_zmem0"].to_numpy(dtype=np.float64)
    conf = dfa[conf_cols].to_numpy(dtype=np.float64)

    part_r = partial_correlation(x, y, conf)
    selection_score = abs(part_r) * (n_analyzable / n_total)

    loo_r, pred_full, pred_resid, actual_resid = loo_diagnostic(x, y, conf)
    is_gap = abs(part_r) - loo_r if (np.isfinite(part_r) and np.isfinite(loo_r)) else float("nan")
    penalty = max(0.0, is_gap) if np.isfinite(is_gap) else float("nan")
    adjusted_score = selection_score - penalty if np.isfinite(penalty) else float("nan")

    per_donor = []
    for row, pf, pr, ar in zip(dfa.itertuples(index=False), pred_full, pred_resid, actual_resid):
        per_donor.append(
            {
                "donor_id": row.donor_id,
                "outcome": float(row.slope_zmem0),
                "predicted": float(pf),
                "feature_value": float(getattr(row, feature_col)),
                "predicted_residual": float(pr),
                "actual_residual": float(ar),
                "abs_error": float(abs(row.slope_zmem0 - pf)),
            }
        )

    return {
        "variation_name": variation_name,
        "feature_name": spec["feature_name"],
        "feature_column": feature_col,
        "partial_r": float(part_r),
        "selection_score": float(selection_score),
        "loo_predictive_r": float(loo_r),
        "is_loo_gap": float(is_gap),
        "penalty": float(penalty),
        "adjusted_score": float(adjusted_score),
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "per_donor_loo": per_donor,
    }


def fmt(x: float) -> str:
    return "nan" if x is None or not np.isfinite(x) else f"{x:.4f}"


def update_canonical_constants(script_path: Path, best_variation: str):
    spec = VARIATION_SPECS[best_variation]
    text = script_path.read_text()
    replacements = {
        "CANONICAL_VARIATION": best_variation,
        "FEATURE_NAME": spec["feature_name"],
        "FEATURE_COLUMN": spec["feature_column"],
    }
    for key, value in replacements.items():
        text = re.sub(rf"^{key}\s*=.*$", f'{key} = "{value}"', text, flags=re.MULTILINE)
    script_path.write_text(text)


def summarize_error_pattern(best_per_donor: list[dict], cohort: pd.DataFrame) -> str:
    if not best_per_donor:
        return "No analyzable donors."
    err = pd.DataFrame(best_per_donor).sort_values("abs_error", ascending=False).head(3)
    extra_cols = [c for c in ["donor_id", "cognitive_status", "braak_numeric", "cerad_ordinal", "sex"] if c in cohort.columns]
    merged = err.merge(cohort[extra_cols], on="donor_id", how="left")
    donors = ", ".join(merged["donor_id"].tolist())
    if "cognitive_status" in merged.columns:
        dementia_n = int((merged["cognitive_status"] == "Dementia").sum())
        dementia_part = f"{dementia_n}/3 are dementia donors and "
    else:
        dementia_part = ""
    if "braak_numeric" in merged.columns:
        high_braak_n = int((merged["braak_numeric"] >= 5).sum())
        braak_part = f"{high_braak_n}/3 have Braak >=5."
    else:
        braak_part = "their shared pathology labels are not all available in the cohort table."
    return f"Largest absolute LOO errors were {donors}; {dementia_part}{braak_part}"


def write_report(
    out_path: Path,
    cohort: pd.DataFrame,
    ranked: list[dict],
    best: dict,
):
    winner = best["variation_name"]
    loser_lines = []
    for r in ranked:
        if r["variation_name"] == winner:
            continue
        loser_lines.append(
            f"{r['variation_name']} (selection {fmt(r['selection_score'])}, partial r {fmt(r['partial_r'])}, LOO {fmt(r['loo_predictive_r'])})"
        )
    loser_summary = "; ".join(loser_lines) if loser_lines else "No alternate local variants were tested."

    if winner == "candidate_variant_b":
        failed_reason = (
            "The majority-reactive threshold looked too permissive: it counts mixed cuffs, "
            "so it likely dilutes the neuron-centered phenotype that stronger reactive dominance is isolating."
        )
        bio_why = (
            "The stricter supermajority threshold better isolates pyramidal neurons whose full local astrocyte cuff "
            "has shifted toward a reactive phenotype."
        )
        observable = "CA1 pyramidal neurons wrapped by astrocyte cuffs in which reactive astrocytes clearly outnumber quiescent astrocytes."
    else:
        failed_reason = (
            "The supermajority threshold looked too strict: it likely makes the event too sparse, "
            "discarding mixed but still meaningfully reactive cuffs."
        )
        bio_why = (
            "The majority-reactive threshold retains enough cuffed neurons to measure a stable donor-level burden "
            "while still requiring the local astrocyte neighborhood to tilt reactive."
        )
        observable = "CA1 pyramidal neurons frequently bordered by astrocyte cuffs where reactive astrocytes are the local majority."

    error_pattern = summarize_error_pattern(best["per_donor_loo"], cohort)

    top_err = pd.DataFrame(best["per_donor_loo"]).sort_values("abs_error", ascending=False).head(3)
    if len(top_err):
        err_examples = ", ".join(
            f"{r.donor_id} (obs {r.outcome:.2f}, pred {r.predicted:.2f})" for r in top_err.itertuples(index=False)
        )
    else:
        err_examples = "No analyzable donors."

    lines = [
        "## Summary",
        f"Tested the {HYPOTHESIS_FAMILY} family; {winner} won with selection score {fmt(best['selection_score'])}.",
        "",
        "## Metrics",
        f"Winner {winner}: partial r {fmt(best['partial_r'])}, selection score {fmt(best['selection_score'])}, "
        f"LOO predictive r {fmt(best['loo_predictive_r'])}, IS-LOO gap {fmt(best['is_loo_gap'])}, "
        f"penalty {fmt(best['penalty'])}, adjusted score {fmt(best['adjusted_score'])}.",
        f"Other tested variations: {loser_summary}",
        "",
        "## Findings",
        "1. What worked and why",
        f"   The winning feature measured the fraction of CA1 pyramidal neurons whose 35 um local astrocyte cuff had at least {MIN_TOTAL_ASTRO} astrocytes and met the reactive-dominance rule. {bio_why}",
        "2. What failed and why",
        f"   {failed_reason}",
        "3. Error pattern: which donors are consistently wrong and what they share",
        f"   {error_pattern} Example high-error donors: {err_examples}",
        "",
        "## Rationale",
        "This approach is biologically coherent because it stays within CA1, keeps the neuron as the center of the niche, and changes only the composition rule of the surrounding cuff.",
        "Compared with nearby alternatives, the winner best balanced specificity against sparsity: it sharpened the cuff phenotype without making the donor score too unstable.",
        "Because the current panel already contains broad astrocytic reactivity and count-based cuff features, this reactive-dominance score is most plausibly additive if cuff composition carries information beyond raw reactive-neighbor count; that final judgment belongs to the panel evaluator.",
        "",
        "## Interpretation",
        f"The signal appears to reflect CA1 peripyramidal astrocyte niches in which the local astrocyte population has shifted toward reactive phenotype dominance around pyramidal neurons.",
        f"Population: CA1 pyramidal neurons with neighboring CA1 astrocytes/reactive astrocytes. Niche: 35 um peripyramidal cuff. Feature summary: donor-level fraction of neurons meeting the reactive-dominant cuff rule. Simplest observable pattern: {observable}",
        "",
        "## Next",
        "Run the next local sweep by keeping the same CA1 peripyramidal niche but varying the minimum cuff size above 2 neighbors, since the current comparison mainly tests composition threshold and the remaining error may come from weakly populated cuffs.",
        "",
    ]
    out_path.write_text("\n".join(lines))


def main():
    data_root = Path("/data")
    scratch = Path("/scratch")
    cohort = load_training_cohort(data_root).copy()

    feature_frames = []
    slide_details = []
    for row in cohort.itertuples(index=False):
        slide_path = data_root / row.slide_name
        slide_result = compute_slide_variations(slide_path)
        rec = {"donor_id": row.donor_id, "slide_name": row.slide_name}
        rec.update({VARIATION_SPECS[k]["feature_column"]: v for k, v in slide_result["scores"].items()})
        rec["n_ca1_pyramidal"] = slide_result["n_ca1_pyramidal"]
        rec["n_ca1_astro_total"] = slide_result["n_ca1_astro_total"]
        rec["n_ca1_reactive_astro"] = slide_result["n_ca1_reactive_astro"]
        feature_frames.append(rec)
        slide_details.append({"donor_id": row.donor_id, **slide_result})

    feat_df = pd.DataFrame(feature_frames)
    donor_df = cohort.merge(feat_df, on=["donor_id", "slide_name"], how="left")
    donor_feature_table_path = scratch / "donor_feature_table.csv"
    donor_df.to_csv(donor_feature_table_path, index=False)

    ranked = [analyze_variation(donor_df, name) for name in VARIATION_SPECS]
    ranked.sort(
        key=lambda d: (
            -np.nan_to_num(d["selection_score"], nan=-np.inf),
            -np.nan_to_num(abs(d["partial_r"]), nan=-np.inf),
            -np.nan_to_num(d["loo_predictive_r"], nan=-np.inf),
        )
    )
    best = ranked[0]

    update_canonical_constants(Path(__file__), best["variation_name"])

    results_payload = {
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": best["variation_name"],
        "feature_name": best["feature_name"],
        "feature_column": best["feature_column"],
        "partial_r": best["partial_r"],
        "selection_score": best["selection_score"],
        "loo_predictive_r": best["loo_predictive_r"],
        "is_loo_gap": best["is_loo_gap"],
        "penalty": best["penalty"],
        "adjusted_score": best["adjusted_score"],
        "n_total": best["n_total"],
        "n_analyzable": best["n_analyzable"],
        "ranked_variations": ranked,
        "slide_details": slide_details,
        "donor_feature_table": str(donor_feature_table_path),
    }
    with open(scratch / "results.json", "w") as f:
        json.dump(results_payload, f, indent=2)

    write_report(scratch / "report.md", cohort, ranked, best)

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best['variation_name']}")
    print(f"  IS partial r:      {fmt(best['partial_r'])}")
    print(f"  Selection score:   {fmt(best['selection_score'])}")
    print(f"  LOO predictive r:  {fmt(best['loo_predictive_r'])}  (diagnostic)")
    print(f"  IS-LOO Gap:        {fmt(best['is_loo_gap'])}  (penalty={fmt(best['penalty'])})")
    print(f"  Adjusted Score:    {fmt(best['adjusted_score'])}")
    print("")
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for r in ranked:
        print(
            f"  {r['variation_name']}  {fmt(r['partial_r'])}  {fmt(r['selection_score'])}  {fmt(r['loo_predictive_r'])}"
        )
    print("")
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best['feature_column']}")
    for rec in best["per_donor_loo"]:
        print(
            f"  {rec['donor_id']}  {rec['outcome']:.4f}  {rec['predicted']:.4f}  {rec['feature_value']:.6f}"
        )


if __name__ == "__main__":
    main()
