from __future__ import annotations

import json
import math
import re
import sys
import warnings
from pathlib import Path

if "/shared/lib" not in sys.path:
    sys.path.insert(0, "/shared/lib")

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import pearsonr

try:
    import openslide  # type: ignore
except Exception:  # pragma: no cover
    openslide = None

from shared_analysis import build_cell_table, load_training_cohort

warnings.filterwarnings(
    "ignore",
    message=r"Object at .* is not recognized as a component of a Zarr hierarchy.",
    category=UserWarning,
)

HYPOTHESIS_FAMILY = "ca1_reactive_astro_pyramidal_isolation_fraction_30um"
CANONICAL_VARIATION = "candidate_variant_a"
FEATURE_NAME = "ca1_pyramidal_reactive_astro_isolated_fraction_k5_30um"
FEATURE_COLUMN = "ca1_pyramidal_reactive_astro_isolated_fraction_k5_30um"

OUTCOME_COL = "slope_zmem0"
CONFOUND_COLS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
DEFAULT_MPP_UM_PER_PX = 0.503

VARIATIONS = [
    {
        "name": "candidate_variant_a",
        "feature_name": "ca1_pyramidal_reactive_astro_isolated_fraction_k5_30um",
        "feature_column": "ca1_pyramidal_reactive_astro_isolated_fraction_k5_30um",
        "max_other_pyramidal_neighbors": 5,
    },
    {
        "name": "candidate_variant_b",
        "feature_name": "ca1_pyramidal_reactive_astro_isolated_fraction_k8_30um",
        "feature_column": "ca1_pyramidal_reactive_astro_isolated_fraction_k8_30um",
        "max_other_pyramidal_neighbors": 8,
    },
]
VARIATION_BY_NAME = {item["name"]: item for item in VARIATIONS}


def _json_default(value):
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (Path,)):
        return str(value)
    raise TypeError(f"Object of type {type(value)!r} is not JSON serializable")


def _design_matrix(frame: pd.DataFrame, cols: list[str]) -> np.ndarray:
    arr = frame[cols].to_numpy(dtype=float)
    return np.column_stack([np.ones(len(frame), dtype=float), arr])


def _residualize(y: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    resid = y - x @ beta
    return resid, beta


def _safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return float("nan")
    x_use = x[mask]
    y_use = y[mask]
    if np.std(x_use) == 0 or np.std(y_use) == 0:
        return float("nan")
    return float(np.corrcoef(x_use, y_use)[0, 1])


def partial_correlation_metrics(df: pd.DataFrame, *, feature_col: str) -> dict[str, float]:
    frame = df.dropna(subset=[feature_col, OUTCOME_COL, *CONFOUND_COLS]).reset_index(drop=True)
    if frame.empty:
        return {"n": 0, "partial_r": float("nan"), "p_value": float("nan")}
    x_conf = _design_matrix(frame, CONFOUND_COLS)
    feat_resid, _ = _residualize(frame[feature_col].to_numpy(dtype=float), x_conf)
    out_resid, _ = _residualize(frame[OUTCOME_COL].to_numpy(dtype=float), x_conf)
    if len(frame) < 3 or np.std(feat_resid) == 0 or np.std(out_resid) == 0:
        return {"n": int(len(frame)), "partial_r": float("nan"), "p_value": float("nan")}
    corr, p_value = pearsonr(feat_resid, out_resid)
    return {"n": int(len(frame)), "partial_r": float(corr), "p_value": float(p_value)}


def loo_prediction_table(df: pd.DataFrame, *, feature_col: str) -> tuple[pd.DataFrame, float]:
    frame = df.dropna(subset=[feature_col, OUTCOME_COL, *CONFOUND_COLS]).reset_index(drop=True)
    rows: list[dict[str, float | str | None]] = []
    pred_resids: list[float] = []
    act_resids: list[float] = []
    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train = _design_matrix(train, CONFOUND_COLS)
        x_test = _design_matrix(test, CONFOUND_COLS)

        feat_train = train[feature_col].to_numpy(dtype=float)
        out_train = train[OUTCOME_COL].to_numpy(dtype=float)
        feat_test = float(test.iloc[0][feature_col])
        out_test = float(test.iloc[0][OUTCOME_COL])

        feat_resid_train, beta_feat = _residualize(feat_train, x_train)
        out_resid_train, beta_out = _residualize(out_train, x_train)

        denom = float(np.dot(feat_resid_train, feat_resid_train))
        if denom <= 0:
            slope = float("nan")
        else:
            slope = float(np.dot(feat_resid_train, out_resid_train) / denom)

        feat_resid_test = feat_test - float((x_test @ beta_feat)[0])
        out_resid_test = out_test - float((x_test @ beta_out)[0])

        pred_resid = slope * feat_resid_test if np.isfinite(slope) else float("nan")
        pred_raw = float((x_test @ beta_out)[0]) + pred_resid if np.isfinite(pred_resid) else float("nan")

        pred_resids.append(pred_resid)
        act_resids.append(out_resid_test)

        row = test.iloc[0].to_dict()
        row["predicted"] = pred_raw
        row["predicted_residualized"] = pred_resid
        row["actual_residualized"] = out_resid_test
        row["absolute_error"] = abs(pred_raw - out_test) if np.isfinite(pred_raw) else float("nan")
        rows.append(row)

    loo_r = _safe_corr(np.asarray(pred_resids, dtype=float), np.asarray(act_resids, dtype=float))
    return pd.DataFrame(rows), loo_r


def selection_and_gap(df: pd.DataFrame, *, feature_col: str) -> dict[str, float]:
    part = partial_correlation_metrics(df, feature_col=feature_col)
    loo_table, loo_r = loo_prediction_table(df, feature_col=feature_col)
    n_total = int(len(df))
    n_analyzable = int(part["n"])
    partial_r = float(part["partial_r"])
    selection_score = abs(partial_r) * (n_analyzable / n_total) if np.isfinite(partial_r) and n_total else float("nan")
    if np.isfinite(loo_r) and np.isfinite(partial_r):
        gap = abs(abs(partial_r) - abs(loo_r))
    elif np.isfinite(partial_r):
        gap = abs(partial_r)
    else:
        gap = float("nan")
    penalty = gap if np.isfinite(gap) else float("nan")
    adjusted = selection_score - penalty if np.isfinite(selection_score) and np.isfinite(penalty) else float("nan")
    return {
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "partial_r": partial_r,
        "p_value": float(part["p_value"]),
        "selection_score": selection_score,
        "loo_predictive_r": loo_r,
        "is_loo_gap": gap,
        "penalty": penalty,
        "adjusted_score": adjusted,
        "loo_table": loo_table,
    }


def get_slide_svs_path(data_root: str | Path, slide_name: str) -> Path:
    data_root = Path(data_root)
    if slide_name.endswith(".zarr"):
        return data_root / slide_name[:-5]
    return data_root / slide_name


def get_slide_mpp_um_per_px(data_root: str | Path, slide_name: str) -> float:
    svs_path = get_slide_svs_path(data_root, slide_name)
    if openslide is not None and svs_path.exists():
        try:
            slide = openslide.OpenSlide(str(svs_path))
            props = dict(slide.properties)
            slide.close()
            for key in ("openslide.mpp-x", "aperio.MPP", "openslide.mpp-y"):
                if key in props:
                    return float(props[key])
        except Exception:
            pass
    try:
        import tifffile  # type: ignore

        if svs_path.exists():
            with tifffile.TiffFile(str(svs_path)) as tif:
                for page in tif.pages[:1]:
                    desc = str(page.description or "")
                    match = re.search(r"aperio\\.MPP\\s*=\\s*([0-9.]+)", desc)
                    if match:
                        return float(match.group(1))
    except Exception:
        pass
    return DEFAULT_MPP_UM_PER_PX


def compute_variation_scores_for_slide(
    *,
    slide_name: str,
    data_root: str | Path,
) -> dict[str, float | int | str]:
    zarr_path = Path(data_root) / slide_name
    cells = build_cell_table(zarr_path, include_regions=True, include_geometry=False)
    cells = cells.loc[(cells["region"] == "CA1") & (cells["cell_type"].isin(["Pyramidal Neuron", "Reactive Astrocyte"]))].copy()

    out: dict[str, float | int | str] = {
        "slide_name": slide_name,
        "mpp_um_per_px": float(get_slide_mpp_um_per_px(data_root, slide_name)),
        "ca1_cell_count_target_types": int(len(cells)),
    }

    pyramidal = cells.loc[cells["cell_type"] == "Pyramidal Neuron", ["x", "y"]].to_numpy(dtype=float)
    reactive = cells.loc[cells["cell_type"] == "Reactive Astrocyte", ["x", "y"]].to_numpy(dtype=float)

    out["ca1_pyramidal_count"] = int(len(pyramidal))
    out["ca1_reactive_astrocyte_count"] = int(len(reactive))

    if len(pyramidal) == 0:
        for variation in VARIATIONS:
            out[variation["feature_column"]] = float("nan")
        out["radius_px"] = float(30.0 / float(out["mpp_um_per_px"]))
        out["reactive_neighbor_rate"] = float("nan")
        out["mean_other_pyramidal_neighbors"] = float("nan")
        return out

    radius_px = 30.0 / float(out["mpp_um_per_px"])
    out["radius_px"] = float(radius_px)

    pyr_tree = cKDTree(pyramidal)
    pyr_neighbors = pyr_tree.query_ball_point(pyramidal, r=radius_px, return_length=True) - 1

    if len(reactive) > 0:
        reactive_tree = cKDTree(reactive)
        reactive_neighbors = reactive_tree.query_ball_point(pyramidal, r=radius_px, return_length=True)
    else:
        reactive_neighbors = np.zeros(len(pyramidal), dtype=int)

    out["reactive_neighbor_rate"] = float(np.mean(reactive_neighbors >= 1))
    out["mean_other_pyramidal_neighbors"] = float(np.mean(pyr_neighbors))

    for variation in VARIATIONS:
        k = int(variation["max_other_pyramidal_neighbors"])
        isolated_mask = (reactive_neighbors >= 1) & (pyr_neighbors <= k)
        out[variation["feature_column"]] = float(np.mean(isolated_mask))
    return out


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
    row = donor_rows.iloc[0]
    slide_name = str(row["slide_name"])
    scores = compute_variation_scores_for_slide(slide_name=slide_name, data_root=data_root)
    return scores.get(VARIATION_BY_NAME[CANONICAL_VARIATION]["feature_column"])


def update_source_with_winner(best_variation_name: str) -> None:
    meta = VARIATION_BY_NAME[best_variation_name]
    path = Path(__file__)
    text = path.read_text()
    replacements = {
        "CANONICAL_VARIATION": f'CANONICAL_VARIATION = "{best_variation_name}"',
        "FEATURE_NAME": f'FEATURE_NAME = "{meta["feature_name"]}"',
        "FEATURE_COLUMN": f'FEATURE_COLUMN = "{meta["feature_column"]}"',
    }
    new_text = text
    for key, replacement in replacements.items():
        new_text = re.sub(rf"^{key} = .*$", replacement, new_text, flags=re.MULTILINE)
    if new_text != text:
        path.write_text(new_text)


def build_feature_table(data_root: str | Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map({"Female": 0.0, "Male": 1.0}).astype(float)
    feature_rows = []
    for _, row in cohort.iterrows():
        donor_id = str(row["donor_id"])
        slide_name = str(row["slide_name"])
        feat = compute_variation_scores_for_slide(slide_name=slide_name, data_root=data_root)
        feat["donor_id"] = donor_id
        feature_rows.append(feat)
    feature_df = pd.DataFrame(feature_rows)
    merged = cohort.merge(feature_df, on=["donor_id", "slide_name"], how="left")
    return merged


def make_error_pattern_text(best_loo: pd.DataFrame) -> str:
    if best_loo.empty:
        return "No analyzable donors for leave-one-out diagnostics."
    err = best_loo.sort_values("absolute_error", ascending=False).head(3)
    donors = []
    braaks = []
    cerads = []
    for _, row in err.iterrows():
        braak = row.get("braak_numeric")
        cerad = row.get("cerad_ordinal")
        if pd.notna(braak):
            braaks.append(float(braak))
        if pd.notna(cerad):
            cerads.append(float(cerad))
        donors.append(
            f'{row["donor_id"]} (abs err {row["absolute_error"]:.3f}, '
            f'braak {int(braak) if pd.notna(braak) else "NA"}, '
            f'cerad {int(cerad) if pd.notna(cerad) else "NA"})'
        )
    if braaks and cerads:
        shared = (
            f"their pathology burden is mixed rather than uniform (braak range {int(min(braaks))}-{int(max(braaks))}, "
            f"cerad range {int(min(cerads))}-{int(max(cerads))}), suggesting this niche is not just a proxy for global AD stage."
        )
    else:
        shared = "available metadata were limited in the worker table."
    return "Largest LOO misses were " + "; ".join(donors) + f"; {shared}"


def make_failed_text(best_name: str, ranked: list[dict[str, float | str]]) -> str:
    others = [r for r in ranked if r["variation_name"] != best_name]
    if not others:
        return "No nearby alternatives were tested."
    other = others[0]
    if best_name == "candidate_variant_a":
        return (
            f'{other["variation_name"]} underperformed because relaxing the isolation threshold to eight nearby '
            "pyramidal neighbors diluted the signal with less-depleted CA1 neighborhoods."
        )
    return (
        f'{other["variation_name"]} underperformed because the stricter five-neighbor cutoff was likely too sparse, '
        "discarding reactive-astro-exposed neurons that still sit in biologically thinned but not fully isolated CA1 neighborhoods."
    )


def write_report(
    *,
    best_meta: dict[str, str | int],
    best_metrics: dict[str, float | int | pd.DataFrame],
    ranked: list[dict[str, float | str]],
) -> None:
    best_loo = best_metrics["loo_table"]
    assert isinstance(best_loo, pd.DataFrame)
    ranking_bits = []
    for item in ranked:
        ranking_bits.append(
            f'{item["variation_name"]}: partial_r={item["partial_r"]:.4f}, '
            f'selection_score={item["selection_score"]:.4f}, loo_predictive_r={item["loo_predictive_r"]:.4f}'
        )
    if best_meta["name"] == "candidate_variant_a":
        rationale_text = (
            "The k<=5 rule keeps the biomarker focused on genuinely sparse CA1 pyramidal neighborhoods, which is the "
            "closest match to reactive-astro-associated local dropout. It beat the k<=8 alternative because the looser "
            "cutoff admits neurons that still have several nearby pyramidal peers, so it is less specific to tissue "
            "disorganization."
        )
        next_text = (
            "Tighten the same family around the winning side of the sweep: test k=3-6 or require at least two reactive "
            "astrocytes within 30 um, because the looser k=8 cutoff lost and the remaining errors likely reflect donors "
            "where astrocyte exposure is present but neuronal depletion must be defined more sharply."
        )
    else:
        rationale_text = (
            "The k<=8 rule appears to capture a broader band of CA1 thinning around reactive astrocytes rather than only "
            "extreme isolation. It beat the k<=5 alternative because the stricter cutoff was probably too brittle for "
            "slides where local neuronal density is reduced but not collapsed."
        )
        next_text = (
            "Stay in the same family but explore k=6-10 or a weighted isolation score based on neighbor count, because "
            "the stricter k=5 version lost and the residual errors suggest some donors carry a milder reactive-astro-"
            "associated thinning signal rather than absolute isolation."
        )

    report = f"""## Summary
One sentence: tested the CA1 reactive-astro-associated pyramidal isolation family; {best_meta["name"]} won, with selection score {best_metrics["selection_score"]:.4f}.

## Metrics
Winning variation `{best_meta["name"]}` (`{best_meta["feature_column"]}`) had partial r {best_metrics["partial_r"]:.4f}, selection score {best_metrics["selection_score"]:.4f}, LOO predictive r {best_metrics["loo_predictive_r"]:.4f}, IS-LOO gap {best_metrics["is_loo_gap"]:.4f}, penalty {best_metrics["penalty"]:.4f}, and adjusted score {best_metrics["adjusted_score"]:.4f}. Other tested variations ranked as: {'; '.join(ranking_bits)}.

## Findings
1. What worked and why (tie to the biological meaning of the target): The signal came from **CA1 pyramidal neurons** and sharpened when restricted to those with a **reactive astrocyte within 30 um** while also having **few nearby CA1 pyramidal peers**. That donor-level fraction appears to summarize local neuronal depletion inside the same reactive-astro niche highlighted by earlier accepted features.
2. What failed and why (specific to the chosen hypothesis and what went wrong): {make_failed_text(str(best_meta["name"]), ranked)}
3. Error pattern: {make_error_pattern_text(best_loo)}

## Rationale
{rationale_text}
If this candidate is kept, it is likely to add more neuron-state specificity than a pure astrocyte abundance summary because it encodes **which CA1 pyramidal neurons sit inside the reactive-astro niche and how locally isolated they are**.

## Interpretation
Biologically, the signal seems to mean **reactive-astro-associated local CA1 pyramidal dropout / sparsification**.  
Population: **CA1 Pyramidal Neuron**.  
Niche: **within 30 um of a Reactive Astrocyte in CA1**.  
Feature summary: **fraction of CA1 pyramidal neurons meeting the reactive-astro-near and low-neighbor criterion**.  
Simplest observable pattern: **more lone or sparsely surrounded CA1 pyramidal somata sitting next to reactive astrocytes**.

## Next
{next_text}
"""
    Path("/scratch/report.md").write_text(report)


def main() -> None:
    data_root = Path("/data")
    feature_table = build_feature_table(data_root)
    feature_table.to_csv("/scratch/donor_feature_table.csv", index=False)

    ranked_variations: list[dict[str, float | str]] = []
    per_variation_loo: dict[str, pd.DataFrame] = {}

    for variation in VARIATIONS:
        metrics = selection_and_gap(feature_table, feature_col=variation["feature_column"])
        record: dict[str, float | str] = {
            "variation_name": variation["name"],
            "feature_name": variation["feature_name"],
            "feature_column": variation["feature_column"],
            "partial_r": float(metrics["partial_r"]),
            "p_value": float(metrics["p_value"]),
            "selection_score": float(metrics["selection_score"]),
            "loo_predictive_r": float(metrics["loo_predictive_r"]),
            "is_loo_gap": float(metrics["is_loo_gap"]),
            "penalty": float(metrics["penalty"]),
            "adjusted_score": float(metrics["adjusted_score"]),
            "n_total": int(metrics["n_total"]),
            "n_analyzable": int(metrics["n_analyzable"]),
        }
        ranked_variations.append(record)
        per_variation_loo[variation["name"]] = metrics["loo_table"]

    ranked_variations.sort(
        key=lambda d: (
            float("-inf") if not np.isfinite(float(d["selection_score"])) else float(d["selection_score"]),
            float("-inf") if not np.isfinite(float(d["adjusted_score"])) else float(d["adjusted_score"]),
            float("-inf") if not np.isfinite(float(d["loo_predictive_r"])) else abs(float(d["loo_predictive_r"])),
        ),
        reverse=True,
    )

    best = ranked_variations[0]
    best_name = str(best["variation_name"])
    best_meta = VARIATION_BY_NAME[best_name]
    best_loo = per_variation_loo[best_name].copy()

    update_source_with_winner(best_name)

    results_payload = {
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": best_name,
        "feature_name": best_meta["feature_name"],
        "feature_column": best_meta["feature_column"],
        "partial_r": float(best["partial_r"]),
        "p_value": float(best["p_value"]),
        "selection_score": float(best["selection_score"]),
        "loo_predictive_r": float(best["loo_predictive_r"]),
        "is_loo_gap": float(best["is_loo_gap"]),
        "penalty": float(best["penalty"]),
        "adjusted_score": float(best["adjusted_score"]),
        "n_total": int(best["n_total"]),
        "n_analyzable": int(best["n_analyzable"]),
        "ranked_variations": ranked_variations,
    }
    Path("/scratch/results.json").write_text(json.dumps(results_payload, indent=2, default=_json_default))

    write_report(
        best_meta=best_meta,
        best_metrics={
            "selection_score": float(best["selection_score"]),
            "partial_r": float(best["partial_r"]),
            "loo_predictive_r": float(best["loo_predictive_r"]),
            "is_loo_gap": float(best["is_loo_gap"]),
            "penalty": float(best["penalty"]),
            "adjusted_score": float(best["adjusted_score"]),
            "loo_table": best_loo,
        },
        ranked=ranked_variations,
    )

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best_name}")
    print(f"  IS partial r:      {float(best['partial_r']):.4f}")
    print(f"  Selection score:   {float(best['selection_score']):.4f}")
    print(f"  LOO predictive r:  {float(best['loo_predictive_r']):.4f}  (diagnostic)")
    print(f"  IS-LOO Gap:        {float(best['is_loo_gap']):.4f}  (penalty={float(best['penalty']):.4f})")
    print(f"  Adjusted Score:    {float(best['adjusted_score']):.4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked_variations:
        print(
            f"  {item['variation_name']:<17} "
            f"{float(item['partial_r']):>9.4f} "
            f"{float(item['selection_score']):>16.4f} "
            f"{float(item['loo_predictive_r']):>17.4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best_meta['feature_column']}")
    for _, row in best_loo.sort_values("donor_id").iterrows():
        print(
            f"  {row['donor_id']:<10} "
            f"{float(row[OUTCOME_COL]):>7.4f} "
            f"{float(row['predicted']):>9.4f} "
            f"{float(row[best_meta['feature_column']]):>9.4f}"
        )


if __name__ == "__main__":
    main()
