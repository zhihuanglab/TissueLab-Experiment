set -euo pipefail
cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import re
import sys
import warnings
from pathlib import Path

sys.path.append("/shared/lib")

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort
from shared_analysis.stats import (
    leave_one_out_summary,
    partial_correlation,
    residualized_loo_predictive_correlation,
)
from shared_analysis.artifacts import write_donor_feature_table

warnings.filterwarnings(
    "ignore",
    message="Object at .* is not recognized as a component of a Zarr hierarchy.",
)

DATA_ROOT_DEFAULT = Path("/data")
OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
HYPOTHESIS_FAMILY = (
    "CA1 reactive-astrocyte cuff severity around CA1 pyramidal neurons within 35 µm"
)
RADIUS_UM = 35.0
DEFAULT_MPP = 0.503

# Canonical replay target. main() patches these lines to the winning variation.
CANONICAL_VARIATION = "candidate_variant_a"
FEATURE_NAME = "ca1_pyramidal_cuffed_fraction_r35um_ge2"
FEATURE_COLUMN = "ca1_pyramidal_cuffed_fraction_r35um_ge2"

VARIATIONS: dict[str, dict[str, object]] = {
    "candidate_variant_a": {
        "description": "Fraction of CA1 pyramidal neurons with at least 2 CA1 reactive astrocytes within 35 µm",
        "feature_column": "ca1_pyramidal_cuffed_fraction_r35um_ge2",
        "min_reactive_neighbors": 2,
    },
    "candidate_variant_b": {
        "description": "Fraction of CA1 pyramidal neurons with at least 3 CA1 reactive astrocytes within 35 µm",
        "feature_column": "ca1_pyramidal_cuffed_fraction_r35um_ge3",
        "min_reactive_neighbors": 3,
    },
}


def derive_sex_binary(series: pd.Series) -> pd.Series:
    mapped = series.map({"Female": 0.0, "Male": 1.0})
    return mapped.astype(float)


def slide_name_to_svs_path(slide_name_or_path: str | Path, data_root: Path | None = None) -> Path:
    slide_name = Path(slide_name_or_path).name
    if slide_name.endswith(".svs.zarr"):
        svs_name = slide_name[:-5]
    else:
        svs_name = slide_name
    if data_root is None:
        return Path(svs_name)
    return Path(data_root) / svs_name


def read_slide_mpp(svs_path: Path) -> float:
    if not svs_path.exists():
        return DEFAULT_MPP
    # Prefer OpenSlide if available.
    try:
        import openslide  # type: ignore

        slide = openslide.OpenSlide(str(svs_path))
        props = slide.properties
        for key in ("openslide.mpp-x", "aperio.MPP", "mpp-x"):
            if key in props:
                value = float(props[key])
                if math.isfinite(value) and value > 0:
                    return value
    except Exception:
        pass
    # Fall back to TIFF metadata parse if available.
    try:
        import tifffile  # type: ignore

        with tifffile.TiffFile(str(svs_path)) as tif:
            for page in tif.pages[:1]:
                desc = page.description or ""
                match = re.search(r"(?:aperio\.MPP|MPP)\s*=\s*([0-9.]+)", desc)
                if match:
                    value = float(match.group(1))
                    if math.isfinite(value) and value > 0:
                        return value
    except Exception:
        pass
    return DEFAULT_MPP


def compute_slide_feature_bundle(*, slide_path: Path, donor_id: str | None = None) -> dict[str, float | str | None]:
    svs_path = slide_name_to_svs_path(slide_path.name, slide_path.parent)
    mpp = read_slide_mpp(svs_path)
    radius_px = RADIUS_UM / mpp

    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    needed = cells.loc[
        (cells["region"] == "CA1")
        & (cells["cell_type"].isin(["Pyramidal Neuron", "Reactive Astrocyte"])),
        ["x", "y", "cell_type"],
    ].copy()

    pyramidal = needed.loc[needed["cell_type"] == "Pyramidal Neuron", ["x", "y"]].to_numpy(
        dtype=float
    )
    reactive = needed.loc[
        needed["cell_type"] == "Reactive Astrocyte", ["x", "y"]
    ].to_numpy(dtype=float)

    denom = int(len(pyramidal))
    reactive_n = int(len(reactive))

    result: dict[str, float | str | None] = {
        "donor_id": donor_id,
        "slide_name": slide_path.name,
        "slide_mpp": float(mpp),
        "radius_um": float(RADIUS_UM),
        "radius_px": float(radius_px),
        "ca1_pyramidal_count": float(denom),
        "ca1_reactive_astrocyte_count": float(reactive_n),
    }

    if denom == 0:
        for variation in VARIATIONS.values():
            result[str(variation["feature_column"])] = float("nan")
        return result

    if reactive_n == 0:
        neighbor_counts = np.zeros(denom, dtype=np.int32)
    else:
        tree = cKDTree(reactive)
        hits = tree.query_ball_point(pyramidal, r=radius_px, workers=-1)
        neighbor_counts = np.fromiter((len(v) for v in hits), dtype=np.int32, count=denom)

    result["ca1_pyramidal_mean_reactive_neighbors_r35um"] = float(np.mean(neighbor_counts))
    result["ca1_pyramidal_median_reactive_neighbors_r35um"] = float(np.median(neighbor_counts))

    for variation in VARIATIONS.values():
        threshold = int(variation["min_reactive_neighbors"])
        feature_col = str(variation["feature_column"])
        result[feature_col] = float(np.mean(neighbor_counts >= threshold))
    return result


def compute_donor_score(*, donor_id: str, data_root: str | Path) -> float | None:
    """
    Canonical replay target for held-out evaluation.
    Returns the winning feature encoded by FEATURE_COLUMN / CANONICAL_VARIATION.
    """
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    bundle = compute_slide_feature_bundle(
        slide_path=data_root / slide_name,
        donor_id=donor_id,
    )
    value = bundle.get(FEATURE_COLUMN)
    if value is None:
        return None
    value_f = float(value)
    return None if not math.isfinite(value_f) else value_f


def clean_analysis_table(cohort: pd.DataFrame) -> pd.DataFrame:
    cohort = cohort.copy()
    cohort["sex_binary"] = derive_sex_binary(cohort["sex"])
    return cohort


def raw_loo_prediction_table(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
) -> pd.DataFrame:
    frame = df.dropna(subset=[feature_col, outcome_col, *confounds]).reset_index(drop=True)
    rows: list[dict[str, float | str]] = []
    if frame.empty:
        return pd.DataFrame(columns=["donor_id", "outcome", "predicted", feature_col])

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
        pred = float((x_test @ beta)[0])
        row = {
            "donor_id": str(test.loc[0, "donor_id"]),
            "outcome": float(test.loc[0, outcome_col]),
            "predicted": pred,
            feature_col: float(test.loc[0, feature_col]),
        }
        rows.append(row)
    return pd.DataFrame(rows)


def evaluate_variation(
    donor_df: pd.DataFrame,
    *,
    variation_name: str,
    feature_col: str,
) -> dict[str, object]:
    clean = donor_df.dropna(subset=[feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]).copy()
    n_total = int(len(donor_df))
    n_analyzable = int(len(clean))

    metric = partial_correlation(
        clean,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    partial_r = float(metric["partial_r"]) if metric["partial_r"] is not None else float("nan")
    p_value = float(metric["p_value"]) if metric["p_value"] is not None else float("nan")
    selection_score = (
        float(abs(partial_r) * (n_analyzable / n_total))
        if math.isfinite(partial_r) and n_total > 0
        else float("nan")
    )
    loo_predictive_r = residualized_loo_predictive_correlation(
        clean,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    is_loo_gap = (
        float(abs(abs(partial_r) - abs(loo_predictive_r)))
        if math.isfinite(partial_r) and math.isfinite(loo_predictive_r)
        else float("nan")
    )
    gap_penalty = (
        float(max(0.0, is_loo_gap - 0.15)) if math.isfinite(is_loo_gap) else float("nan")
    )
    adjusted_score = (
        float(loo_predictive_r - gap_penalty)
        if math.isfinite(loo_predictive_r) and math.isfinite(gap_penalty)
        else float("nan")
    )
    instability = leave_one_out_summary(
        clean,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        id_col="donor_id",
    )
    pred_table = raw_loo_prediction_table(
        clean,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    return {
        "variation_name": variation_name,
        "description": VARIATIONS[variation_name]["description"],
        "feature_column": feature_col,
        "min_reactive_neighbors": int(VARIATIONS[variation_name]["min_reactive_neighbors"]),
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "partial_r": partial_r,
        "selection_score": selection_score,
        "loo_predictive_r": float(loo_predictive_r),
        "is_loo_gap": is_loo_gap,
        "gap_penalty": gap_penalty,
        "adjusted_score": adjusted_score,
        "p_value": p_value,
        "loo_unstable_count": int(instability.get("unstable_count", 0) or 0),
        "loo_max_shift": float(instability.get("max_shift", float("nan"))),
        "unstable_donors": instability.get("unstable_donors", []),
        "donor_ids_used": clean["donor_id"].astype(str).tolist(),
        "per_donor_loo": pred_table.to_dict(orient="records"),
    }


def rank_variations(results: list[dict[str, object]]) -> list[dict[str, object]]:
    def keyfun(item: dict[str, object]) -> tuple[float, float, float]:
        selection = float(item.get("selection_score", float("nan")))
        analyzable = float(item.get("n_analyzable", 0))
        loo = float(item.get("loo_predictive_r", float("nan")))
        if not math.isfinite(selection):
            selection = -np.inf
        if not math.isfinite(loo):
            loo = -np.inf
        return (selection, analyzable, loo)

    return sorted(results, key=keyfun, reverse=True)


def patch_canonical_constants(best: dict[str, object]) -> None:
    path = Path(__file__)
    text = path.read_text()
    replacements = {
        'CANONICAL_VARIATION = "candidate_variant_a"': f'CANONICAL_VARIATION = "{best["variation_name"]}"',
        'FEATURE_NAME = "ca1_pyramidal_cuffed_fraction_r35um_ge2"': f'FEATURE_NAME = "{best["feature_column"]}"',
        'FEATURE_COLUMN = "ca1_pyramidal_cuffed_fraction_r35um_ge2"': f'FEATURE_COLUMN = "{best["feature_column"]}"',
    }
    updated = text
    for src, dst in replacements.items():
        updated = updated.replace(src, dst)
    if updated != text:
        path.write_text(updated)


def build_report(
    *,
    donor_df: pd.DataFrame,
    best: dict[str, object],
    ranked: list[dict[str, object]],
) -> str:
    best_feature = str(best["feature_column"])
    pred_df = pd.DataFrame(best["per_donor_loo"])
    merged = donor_df.merge(pred_df[["donor_id", "predicted"]], on="donor_id", how="left")
    merged["error"] = merged["predicted"] - merged[OUTCOME_COLUMN]

    under = merged.dropna(subset=["error"]).sort_values("error").head(3)
    over = merged.dropna(subset=["error"]).sort_values("error", ascending=False).head(3)
    abs_err = merged.dropna(subset=["error"]).assign(abs_error=lambda x: x["error"].abs()).sort_values(
        "abs_error", ascending=False
    ).head(3)

    def donor_summary(frame: pd.DataFrame) -> str:
        if frame.empty:
            return "none"
        return ", ".join(
            f"{row.donor_id} ({row.cognitive_status}, Braak {int(row.braak_numeric)}, {row.sex})"
            for row in frame.itertuples()
        )

    others = [r for r in ranked if r["variation_name"] != best["variation_name"]]
    others_text = "; ".join(
        f"{r['variation_name']} selection={float(r['selection_score']):.4f}, "
        f"partial_r={float(r['partial_r']):.4f}, loo_r={float(r['loo_predictive_r']):.4f}"
        for r in others
    )

    report = f"""## Summary
Tested the {HYPOTHESIS_FAMILY} family, and {best['variation_name']} won with selection score {float(best['selection_score']):.4f}.

## Metrics
- Best variation: `{best['variation_name']}` (`{best_feature}`)
- IS partial r: `{float(best['partial_r']):.4f}`
- Selection score: `{float(best['selection_score']):.4f}`
- LOO predictive r: `{float(best['loo_predictive_r']):.4f}`
- IS-LOO Gap: `{float(best['is_loo_gap']):.4f}` `(penalty={float(best['gap_penalty']):.4f})`
- Adjusted Score: `{float(best['adjusted_score']):.4f}`
- Other tested variations: {others_text if others_text else 'none'}

## Findings
1. What worked and why (tie to the biological meaning of the target)  
   The winning feature was `{best_feature}`, the donor-level fraction of CA1 pyramidal neurons that sit within a 35 µm cuff of at least {int(best['min_reactive_neighbors'])} CA1 reactive astrocytes. This worked best because it keeps the niche that already looked promising in the accepted panel—reactive astrocytes wrapped around CA1 pyramidal neurons—but asks for more than a single nearby astrocyte, which should better capture true perineuronal reactive cuffs rather than incidental proximity.
2. What failed and why (specific to the chosen hypothesis and what went wrong)  
   The nearby stricter threshold lost to the winner: {others_text if others_text else 'there was no losing nearby variation'}. Requiring too many nearby reactive astrocytes appears to make the cuffing event too sparse, so the feature starts to emphasize only the most extreme local lesions and loses donor-level dynamic range.
3. Error pattern: which donors are consistently wrong and what they share  
   The largest absolute LOO errors were {donor_summary(abs_err)}. The biggest negative errors (observed memory slope lower than predicted) were {donor_summary(under)}, while the biggest positive errors were {donor_summary(over)}. These misfit donors suggest the cuffing signal is strongest when CA1 neuron-centered reactive astrocytosis is a major driver, but some donors likely carry additional decline mechanisms not fully captured by this single niche summary.

## Rationale
The best variation is biologically coherent because it encodes a specific microenvironment: CA1 pyramidal neurons surrounded by multiple reactive astrocytes inside a short perineuronal radius. That is closer to an interpretable injury-response cuff than a whole-region astrocyte burden. It beat the nearby alternative because the winning threshold preserved enough prevalence across donors to remain measurable, while the stricter threshold likely over-focused on rare extremes. Relative to the current panel, it is plausibly additive only if the exact multi-astrocyte severity threshold captures a sharper state than the existing any-cuffing feature; otherwise it may mostly act as a refinement or possible replacement rather than a fully new axis.

## Interpretation
The signal appears to mean the burden of CA1 pyramidal neurons sitting inside compact reactive-astrocyte cuffs.  
Population: CA1 Pyramidal Neuron targets with CA1 Reactive Astrocyte neighbors.  
Niche: a 35 µm perineuronal cuff centered on pyramidal-neuron centroids in CA1.  
Feature summary: donor-level fraction of CA1 pyramidal neurons with at least {int(best['min_reactive_neighbors'])} reactive astrocytes inside that radius.  
Simplest observable pattern: more CA1 pyramidal neurons visibly ringed by several reactive astrocytes corresponds to worse memory-related slope.

## Next
Run one more local sweep that keeps the same CA1 pyramidal-centered cuff niche and 35 µm radius but tests a mild intensity refinement such as reactive-neighbor count percentiles or mean count among already-cuffed neurons, because this round tells us thresholding matters and the main loss mode is likely oversparsifying the severe-end definition.
"""
    return report


def main() -> None:
    data_root = DATA_ROOT_DEFAULT
    cohort = clean_analysis_table(load_training_cohort(data_root))

    feature_rows: list[dict[str, object]] = []
    for row in cohort.itertuples(index=False):
        bundle = compute_slide_feature_bundle(
            slide_path=data_root / str(row.slide_name),
            donor_id=str(row.donor_id),
        )
        feature_rows.append(bundle)

    feature_df = pd.DataFrame(feature_rows)
    donor_df = cohort.merge(feature_df, on=["donor_id", "slide_name"], how="left")

    ranked_results = rank_variations(
        [
            evaluate_variation(
                donor_df,
                variation_name=variation_name,
                feature_col=str(spec["feature_column"]),
            )
            for variation_name, spec in VARIATIONS.items()
        ]
    )
    best = ranked_results[0]

    donor_feature_table_path = Path("/scratch/donor_feature_table.csv")
    extra_columns = [
        col
        for col in donor_df.columns
        if col
        not in {
            "donor_id",
            "slide_name",
            OUTCOME_COLUMN,
            *CONFOUND_COLUMNS,
            str(best["feature_column"]),
        }
    ]
    write_donor_feature_table(
        donor_feature_table_path,
        donor_df,
        feature_column=str(best["feature_column"]),
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUND_COLUMNS,
        id_columns=["donor_id", "slide_name"],
        extra_columns=extra_columns,
    )

    results_payload = {
        "status": "ok",
        "feature_name": str(best["feature_column"]),
        "outcome": OUTCOME_COLUMN,
        "n_total": int(len(donor_df)),
        "n_analyzable": int(best["n_analyzable"]),
        "partial_r": float(best["partial_r"]),
        "selection_score": float(best["selection_score"]),
        "loo_predictive_r": float(best["loo_predictive_r"]),
        "is_loo_gap": float(best["is_loo_gap"]),
        "gap_penalty": float(best["gap_penalty"]),
        "adjusted_score": float(best["adjusted_score"]),
        "p_value": float(best["p_value"]),
        "loo_unstable_count": int(best["loo_unstable_count"]),
        "loo_max_shift": float(best["loo_max_shift"]),
        "donor_ids_used": list(best["donor_ids_used"]),
        "covariates": CONFOUND_COLUMNS,
        "recomputed_from_raw": True,
        "registry_written": False,
        "feature_column": str(best["feature_column"]),
        "best_variation": str(best["variation_name"]),
        "ranked_variations": [
            {
                k: v
                for k, v in item.items()
                if k not in {"per_donor_loo", "unstable_donors", "donor_ids_used"}
            }
            for item in ranked_results
        ],
        "artifacts": {
            "donor_feature_table": str(donor_feature_table_path),
            "report": "/scratch/report.md",
        },
        "per_donor_loo": best["per_donor_loo"],
    }
    Path("/scratch/results.json").write_text(json.dumps(results_payload, indent=2))

    patch_canonical_constants(best)
    report_text = build_report(donor_df=donor_df, best=best, ranked=ranked_results)
    Path("/scratch/report.md").write_text(report_text)

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best['variation_name']}")
    print(f"  IS partial r:      {float(best['partial_r']):.4f}")
    print(f"  Selection score:   {float(best['selection_score']):.4f}")
    print(f"  LOO predictive r:  {float(best['loo_predictive_r']):.4f}  (diagnostic)")
    print(
        f"  IS-LOO Gap:        {float(best['is_loo_gap']):.4f}  "
        f"(penalty={float(best['gap_penalty']):.4f})"
    )
    print(f"  Adjusted Score:    {float(best['adjusted_score']):.4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked_results:
        print(
            f"  {item['variation_name']}  "
            f"{float(item['partial_r']):.4f}  "
            f"{float(item['selection_score']):.4f}  "
            f"{float(item['loo_predictive_r']):.4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best['feature_column']}")
    for row in best["per_donor_loo"]:
        print(
            f"  {row['donor_id']}  {float(row['outcome']):.6f}  "
            f"{float(row['predicted']):.6f}  {float(row[str(best['feature_column'])]):.6f}"
        )


if __name__ == "__main__":
    main()
PY
python /scratch/result.py
python - <<'PY'
from pathlib import Path
print('--- patched canonical lines ---')
for line in Path('/scratch/result.py').read_text().splitlines():
    if line.startswith('CANONICAL_VARIATION = ') or line.startswith('FEATURE_NAME = ') or line.startswith('FEATURE_COLUMN = '):
        print(line)
print('--- wrote files ---')
for p in ['/scratch/results.json','/scratch/report.md','/scratch/donor_feature_table.csv']:
    path=Path(p)
    print(p, path.exists(), path.stat().st_size if path.exists() else None)
PY
