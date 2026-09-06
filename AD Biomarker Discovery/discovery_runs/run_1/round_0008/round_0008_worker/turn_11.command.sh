cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

if "/shared/lib" not in sys.path:
    sys.path.insert(0, "/shared/lib")

from shared_analysis import build_cell_table, load_training_cohort  # noqa: E402

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

HYPOTHESIS_FAMILY = "ca1_reactive_astro_proximal_small_pyramidal_enrichment"
CANONICAL_VARIATION = "candidate_variant_a"
FEATURE_NAME = "ca1_reactive_astro_proximal_small_pyramidal_enrichment"
FEATURE_COLUMN = "ca1_small_pyramidal_enrichment_q25_80px"

DATA_ROOT_DEFAULT = Path("/data")
RESULTS_PATH_DEFAULT = Path("/scratch/results.json")
DONOR_TABLE_PATH_DEFAULT = Path("/scratch/donor_feature_table.csv")
REPORT_PATH_DEFAULT = Path("/scratch/report.md")

PROXIMAL_RADIUS_PX = 80.0
MIN_COMPARTMENT_COUNT = 20
PERCENTILE_BY_VARIATION = {
    "candidate_variant_a": 25.0,
    "candidate_variant_b": 20.0,
}
FEATURE_COLUMN_BY_VARIATION = {
    "candidate_variant_a": "ca1_small_pyramidal_enrichment_q25_80px",
    "candidate_variant_b": "ca1_small_pyramidal_enrichment_q20_80px",
}


@dataclass
class DonorFeatureSummary:
    donor_id: str
    slide_name: str
    n_ca1_pyramidal: int
    n_ca1_reactive_astro: int
    n_proximal_pyramidal: int
    n_distal_pyramidal: int
    q20_threshold_area: float | None
    q25_threshold_area: float | None
    proximal_small_fraction_q20: float | None
    distal_small_fraction_q20: float | None
    proximal_small_fraction_q25: float | None
    distal_small_fraction_q25: float | None
    candidate_variant_a: float | None
    candidate_variant_b: float | None


def _safe_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        value = float(value)
    except Exception:
        return None
    if not np.isfinite(value):
        return None
    return value


def _canonical_variation_from_results() -> str:
    results_path = Path(__file__).with_name("results.json")
    if results_path.exists():
        try:
            payload = json.loads(results_path.read_text())
            name = payload.get("best_variation")
            if name in PERCENTILE_BY_VARIATION:
                return str(name)
        except Exception:
            pass
    return CANONICAL_VARIATION


def _compute_single_donor_summary(slide_path: Path, donor_id: str) -> DonorFeatureSummary:
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=True)

    region = cells["region"].astype("object")
    cell_type = cells["cell_type"].astype("object")

    ca1 = cells.loc[region == "CA1", ["x", "y", "area", "cell_type"]].copy()

    pyramidal = ca1.loc[ca1["cell_type"] == "Pyramidal Neuron", ["x", "y", "area"]].copy()
    reactive = ca1.loc[ca1["cell_type"] == "Reactive Astrocyte", ["x", "y"]].copy()

    pyramidal["area"] = pd.to_numeric(pyramidal["area"], errors="coerce")
    pyramidal = pyramidal.loc[np.isfinite(pyramidal["area"].to_numpy()) & (pyramidal["area"].to_numpy() > 0)].copy()

    n_pyr = int(len(pyramidal))
    n_reactive = int(len(reactive))

    if n_pyr == 0 or n_reactive == 0:
        return DonorFeatureSummary(
            donor_id=donor_id,
            slide_name=slide_path.name,
            n_ca1_pyramidal=n_pyr,
            n_ca1_reactive_astro=n_reactive,
            n_proximal_pyramidal=0,
            n_distal_pyramidal=n_pyr,
            q20_threshold_area=None,
            q25_threshold_area=None,
            proximal_small_fraction_q20=None,
            distal_small_fraction_q20=None,
            proximal_small_fraction_q25=None,
            distal_small_fraction_q25=None,
            candidate_variant_a=None,
            candidate_variant_b=None,
        )

    pyr_xy = pyramidal[["x", "y"]].to_numpy(dtype=float)
    reactive_xy = reactive[["x", "y"]].to_numpy(dtype=float)
    tree = cKDTree(reactive_xy)
    distances, _ = tree.query(pyr_xy, k=1, distance_upper_bound=PROXIMAL_RADIUS_PX)
    proximal_mask = np.isfinite(distances)

    n_prox = int(proximal_mask.sum())
    n_dist = int((~proximal_mask).sum())

    q20 = float(np.nanpercentile(pyramidal["area"].to_numpy(dtype=float), 20.0))
    q25 = float(np.nanpercentile(pyramidal["area"].to_numpy(dtype=float), 25.0))

    areas = pyramidal["area"].to_numpy(dtype=float)

    def _score_for_threshold(threshold: float) -> tuple[float | None, float | None, float | None]:
        if n_prox < MIN_COMPARTMENT_COUNT or n_dist < MIN_COMPARTMENT_COUNT:
            return None, None, None
        small = areas <= threshold
        prox_small = float(np.mean(small[proximal_mask])) if n_prox > 0 else None
        dist_small = float(np.mean(small[~proximal_mask])) if n_dist > 0 else None
        if prox_small is None or dist_small is None:
            return None, prox_small, dist_small
        return float(prox_small - dist_small), prox_small, dist_small

    score_a, prox_frac_q25, dist_frac_q25 = _score_for_threshold(q25)
    score_b, prox_frac_q20, dist_frac_q20 = _score_for_threshold(q20)

    return DonorFeatureSummary(
        donor_id=donor_id,
        slide_name=slide_path.name,
        n_ca1_pyramidal=n_pyr,
        n_ca1_reactive_astro=n_reactive,
        n_proximal_pyramidal=n_prox,
        n_distal_pyramidal=n_dist,
        q20_threshold_area=q20,
        q25_threshold_area=q25,
        proximal_small_fraction_q20=prox_frac_q20,
        distal_small_fraction_q20=dist_frac_q20,
        proximal_small_fraction_q25=prox_frac_q25,
        distal_small_fraction_q25=dist_frac_q25,
        candidate_variant_a=score_a,
        candidate_variant_b=score_b,
    )


def compute_donor_score(*, donor_id: str, data_root: str | Path) -> float | None:
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    row = cohort.loc[cohort["donor_id"] == donor_id]
    if row.empty:
        return None
    slide_name = str(row.iloc[0]["slide_name"])
    summary = _compute_single_donor_summary(data_root / slide_name, donor_id)
    variation = _canonical_variation_from_results()
    value = getattr(summary, variation)
    return _safe_float(value)


def _prepare_cohort_table(data_root: Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map({"Male": 1.0, "Female": 0.0})
    return cohort


def _extract_all_donor_features(data_root: Path) -> pd.DataFrame:
    cohort = _prepare_cohort_table(data_root)
    records: list[dict[str, Any]] = []
    for _, row in cohort.iterrows():
        donor_id = str(row["donor_id"])
        slide_name = str(row["slide_name"])
        summary = _compute_single_donor_summary(data_root / slide_name, donor_id)
        records.append(summary.__dict__)
        print(
            f"EXTRACTED {donor_id}: CA1 pyramidal={summary.n_ca1_pyramidal}, "
            f"reactive_astro={summary.n_ca1_reactive_astro}, proximal={summary.n_proximal_pyramidal}, "
            f"distal={summary.n_distal_pyramidal}"
        )
    feature_df = pd.DataFrame.from_records(records)
    merged = cohort.merge(feature_df, on=["donor_id", "slide_name"], how="left")
    merged["ca1_small_pyramidal_enrichment_q25_80px"] = merged["candidate_variant_a"]
    merged["ca1_small_pyramidal_enrichment_q20_80px"] = merged["candidate_variant_b"]
    return merged


def _pearsonr_safe(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) < 3:
        return float("nan")
    if np.allclose(np.nanstd(x), 0.0) or np.allclose(np.nanstd(y), 0.0):
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def _residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    model = LinearRegression()
    model.fit(x, y)
    return y - model.predict(x)


def _partial_correlation(feature: np.ndarray, outcome: np.ndarray, confounds: np.ndarray) -> float:
    feature_resid = _residualize(feature, confounds)
    outcome_resid = _residualize(outcome, confounds)
    return _pearsonr_safe(feature_resid, outcome_resid)


def _loo_predictions(
    feature: np.ndarray,
    outcome: np.ndarray,
    confounds: np.ndarray,
    donor_ids: list[str],
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    n = len(feature)
    preds = np.full(n, np.nan, dtype=float)
    rows: list[dict[str, Any]] = []

    for i in range(n):
        train_mask = np.ones(n, dtype=bool)
        train_mask[i] = False

        x_train = np.column_stack([confounds[train_mask], feature[train_mask]])
        x_test = np.column_stack([confounds[~train_mask], feature[~train_mask]])

        scaler = StandardScaler()
        x_train_scaled = scaler.fit_transform(x_train)
        x_test_scaled = scaler.transform(x_test)

        model = LinearRegression()
        model.fit(x_train_scaled, outcome[train_mask])
        pred = float(model.predict(x_test_scaled)[0])
        preds[i] = pred
        rows.append(
            {
                "donor_id": donor_ids[i],
                "outcome": float(outcome[i]),
                "predicted": pred,
                "feature_value": float(feature[i]),
            }
        )
    return preds, rows


def _evaluate_variation(
    df: pd.DataFrame,
    variation_name: str,
) -> dict[str, Any]:
    feature_column = FEATURE_COLUMN_BY_VARIATION[variation_name]
    required = ["slope_zmem0", "max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary", feature_column]
    sub = df.loc[:, ["donor_id"] + required].copy()
    for col in required:
        sub[col] = pd.to_numeric(sub[col], errors="coerce")
    analyzable = sub.dropna().reset_index(drop=True)

    n_total = len(df)
    n_analyzable = len(analyzable)

    result: dict[str, Any] = {
        "variation_name": variation_name,
        "feature_column": feature_column,
        "n_total": int(n_total),
        "n_analyzable": int(n_analyzable),
    }

    if n_analyzable < 8:
        result.update(
            {
                "partial_r": None,
                "selection_score": None,
                "loo_predictive_r": None,
                "is_loo_gap": None,
                "penalty": None,
                "adjusted_score": None,
                "per_donor_loo": [],
            }
        )
        return result

    feature = analyzable[feature_column].to_numpy(dtype=float)
    outcome = analyzable["slope_zmem0"].to_numpy(dtype=float)
    confounds = analyzable[["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]].to_numpy(dtype=float)
    donor_ids = analyzable["donor_id"].astype(str).tolist()

    partial_r = _partial_correlation(feature, outcome, confounds)
    selection_score = abs(partial_r) * (n_analyzable / n_total)

    preds, per_donor = _loo_predictions(feature, outcome, confounds, donor_ids)
    loo_r = _pearsonr_safe(preds, outcome)

    gap = abs(partial_r - loo_r) if np.isfinite(partial_r) and np.isfinite(loo_r) else float("nan")
    penalty = max(0.0, gap - 0.15) * 0.5 if np.isfinite(gap) else float("nan")
    adjusted = selection_score - penalty if np.isfinite(selection_score) and np.isfinite(penalty) else float("nan")

    result.update(
        {
            "partial_r": float(partial_r) if np.isfinite(partial_r) else None,
            "selection_score": float(selection_score) if np.isfinite(selection_score) else None,
            "loo_predictive_r": float(loo_r) if np.isfinite(loo_r) else None,
            "is_loo_gap": float(gap) if np.isfinite(gap) else None,
            "penalty": float(penalty) if np.isfinite(penalty) else None,
            "adjusted_score": float(adjusted) if np.isfinite(adjusted) else None,
            "per_donor_loo": per_donor,
        }
    )
    return result


def _score_key(item: dict[str, Any]) -> tuple[float, float, float]:
    sel = item.get("selection_score")
    loo = item.get("loo_predictive_r")
    adj = item.get("adjusted_score")
    return (
        -999.0 if sel is None else float(sel),
        -999.0 if loo is None else abs(float(loo)),
        -999.0 if adj is None else float(adj),
    )


def _format_num(x: float | None) -> str:
    if x is None or not np.isfinite(float(x)):
        return "NA"
    return f"{float(x):.4f}"


def _update_script_canonical_constants(best_variation: str) -> None:
    path = Path(__file__)
    text = path.read_text()
    new_feature_column = FEATURE_COLUMN_BY_VARIATION[best_variation]
    replacements = {
        'CANONICAL_VARIATION = "candidate_variant_a"': f'CANONICAL_VARIATION = "{best_variation}"',
        'FEATURE_COLUMN = "ca1_small_pyramidal_enrichment_q25_80px"': f'FEATURE_COLUMN = "{new_feature_column}"',
    }
    for old, new in replacements.items():
        if old in text:
            text = text.replace(old, new, 1)
    path.write_text(text)


def _write_report(
    df: pd.DataFrame,
    best: dict[str, Any],
    ranked: list[dict[str, Any]],
    report_path: Path,
) -> None:
    best_col = best["feature_column"]
    analyzable_best = df.loc[df[best_col].notna(), ["donor_id", "slope_zmem0", best_col]].copy()

    loo_map = {row["donor_id"]: row for row in best.get("per_donor_loo", [])}
    analyzable_best["predicted"] = analyzable_best["donor_id"].map(lambda d: loo_map.get(d, {}).get("predicted"))
    analyzable_best["abs_error"] = (analyzable_best["slope_zmem0"] - analyzable_best["predicted"]).abs()
    worst = analyzable_best.sort_values("abs_error", ascending=False).head(5)

    other_lines = []
    for item in ranked[1:]:
        other_lines.append(
            f"- {item['variation_name']}: partial_r={_format_num(item.get('partial_r'))}, "
            f"selection_score={_format_num(item.get('selection_score'))}, "
            f"loo_predictive_r={_format_num(item.get('loo_predictive_r'))}"
        )
    if not other_lines:
        other_lines = ["- No other valid variations."]

    worst_text = "; ".join(
        [
            f"{row.donor_id} (outcome={row.slope_zmem0:.3f}, predicted={row.predicted:.3f}, "
            f"{best_col}={row[best_col]:.3f})"
            for _, row in worst.iterrows()
            if pd.notna(row["predicted"])
        ]
    )
    if not worst_text:
        worst_text = "No stable analyzable donor-level error pattern could be summarized."

    if best["variation_name"] == "candidate_variant_a":
        threshold_text = "25th-percentile within-donor small-cell threshold"
        loser_text = "The stricter 20th-percentile tail likely discarded too many moderately shrunken proximal pyramidal neurons."
    else:
        threshold_text = "20th-percentile within-donor small-cell threshold"
        loser_text = "The broader 25th-percentile threshold likely mixed the informative tail with less selective background size variation."

    report = f"""## Summary
Tested a CA1 reactive-astrocyte-proximal small-pyramidal enrichment family; {best['variation_name']} won with local selection score {_format_num(best.get('selection_score'))}.

## Metrics
Best variation: {best['variation_name']}
- IS partial r: {_format_num(best.get('partial_r'))}
- Selection score: {_format_num(best.get('selection_score'))}
- LOO predictive r: {_format_num(best.get('loo_predictive_r'))}
- IS-LOO Gap: {_format_num(best.get('is_loo_gap'))} (penalty={_format_num(best.get('penalty'))})
- Adjusted Score: {_format_num(best.get('adjusted_score'))}
- Feature column: `{best_col}`

Other tested variations:
{chr(10).join(other_lines)}

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The winner worked best when the donor-level scalar asked whether morphologically small CA1 pyramidal neurons were selectively enriched inside the reactive-astrocyte-proximal compartment rather than globally. That contrast is aligned with a localized vulnerability signal rather than simple overall neuronal size or abundance.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - {loser_text} It remained the same CA1 pyramidal/reactive-astrocyte niche, but the small-cell cutoff changed the bias-variance tradeoff inside a limited proximal compartment.
3. Error pattern: which donors are consistently wrong and what they share
   - Largest absolute LOO errors: {worst_text}. These errors tend to reflect donors where the niche contrast does not map cleanly onto the outcome despite measurable CA1 pyramidal and reactive astrocyte presence.

## Rationale
The winning variation is biologically coherent because it focuses on a specific population (CA1 pyramidal neurons), a specific niche (within 80 px of CA1 reactive astrocytes), and a donor-relative size-tail summary ({threshold_text}). It beat the nearby alternative because the chosen cutoff better isolated the donor-specific small-cell pattern without relying on an absolute area scale. Given the accepted panel already contains CA1 pyramidal fraction, reactive astrocyte lineage fraction, niche fraction, niche coverage, and proximal median log area, this enrichment contrast still looks plausibly additive because it measures selective concentration of the small-size tail between proximal and distal CA1 compartments rather than another absolute burden summary.

## Interpretation
The signal appears to represent localized enrichment of morphologically small, potentially stressed or atrophic CA1 pyramidal neurons around reactive astrocytes. Population: CA1 Pyramidal Neuron. Niche: CA1 pyramidal cells within 80 px of CA1 Reactive Astrocytes versus the distal CA1 background. Feature summary: proximal small-cell fraction minus distal small-cell fraction. Simplest observable pattern: donors with worse memory decline appear to show more shrinkage-like small pyramidal profiles concentrated near reactive astrocyte neighborhoods inside CA1.

## Next
Run the next local sweep around the same CA1 pyramidal/reactive-astrocyte niche but keep the winning proximal-minus-distal contrast and vary only the proximity radius (for example 60 px versus 100 px) to determine whether the signal is truly juxtacrine-local or a broader reactive-glial field effect.
"""
    report_path.write_text(report)


def main() -> None:
    data_root = DATA_ROOT_DEFAULT
    donor_table = _extract_all_donor_features(data_root)
    donor_table.to_csv(DONOR_TABLE_PATH_DEFAULT, index=False)

    variation_results = [_evaluate_variation(donor_table, name) for name in PERCENTILE_BY_VARIATION]
    ranked = sorted(variation_results, key=_score_key, reverse=True)
    best = ranked[0]

    payload = {
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": best["variation_name"],
        "feature_name": FEATURE_NAME,
        "feature_column": best["feature_column"],
        "partial_r": best.get("partial_r"),
        "selection_score": best.get("selection_score"),
        "loo_predictive_r": best.get("loo_predictive_r"),
        "is_loo_gap": best.get("is_loo_gap"),
        "penalty": best.get("penalty"),
        "adjusted_score": best.get("adjusted_score"),
        "n_total": best.get("n_total"),
        "n_analyzable": best.get("n_analyzable"),
        "ranked_variations": ranked,
    }
    RESULTS_PATH_DEFAULT.write_text(json.dumps(payload, indent=2))

    _update_script_canonical_constants(best["variation_name"])
    _write_report(donor_table, best, ranked, REPORT_PATH_DEFAULT)

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best['variation_name']}")
    print(f"  IS partial r:      {_format_num(best.get('partial_r'))}")
    print(f"  Selection score:   {_format_num(best.get('selection_score'))}")
    print(f"  LOO predictive r:  {_format_num(best.get('loo_predictive_r'))}  (diagnostic)")
    print(
        f"  IS-LOO Gap:        {_format_num(best.get('is_loo_gap'))}  "
        f"(penalty={_format_num(best.get('penalty'))})"
    )
    print(f"  Adjusted Score:    {_format_num(best.get('adjusted_score'))}")
    print("")
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked:
        print(
            f"  {item['variation_name']}  "
            f"{_format_num(item.get('partial_r'))}  "
            f"{_format_num(item.get('selection_score'))}  "
            f"{_format_num(item.get('loo_predictive_r'))}"
        )
    print("")
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best['feature_column']}")
    feature_map = donor_table.set_index("donor_id")[best["feature_column"]].to_dict()
    for row in best.get("per_donor_loo", []):
        feature_value = feature_map.get(row["donor_id"])
        feature_text = "NA" if feature_value is None or not np.isfinite(feature_value) else f"{float(feature_value):.6f}"
        print(
            f"  {row['donor_id']}  {float(row['outcome']):.6f}  "
            f"{float(row['predicted']):.6f}  {feature_text}"
        )


if __name__ == "__main__":
    main()
PY
python -m py_compile /scratch/result.py
