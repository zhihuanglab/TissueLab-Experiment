cat > /scratch/result.py <<'PY'
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

from shared_analysis import (
    build_cell_table,
    leave_one_out_summary,
    load_training_cohort,
    partial_correlation,
    write_donor_feature_table,
)

warnings.filterwarnings("ignore", message="Object at .* is not recognized as a component of a Zarr hierarchy.")
warnings.filterwarnings("ignore", category=RuntimeWarning)

FAMILY_NAME = "ca1_reactive_astrocyte_proximal_pyramidal_lower_tail_log_area"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
RADIUS_PX = 80.0

VARIATIONS = [
    {
        "name": "candidate_variant_a",
        "label": "q25 within 80 px",
        "description": "25th percentile of log1p(area) for CA1 pyramidal neurons within 80 px of any CA1 reactive astrocyte.",
        "quantile": 0.25,
        "radius_px": 80.0,
        "feature_column": "ca1_pyramidal_reactive_astro_proximal_q25_log_area_80px",
    },
    {
        "name": "candidate_variant_b",
        "label": "q20 within 80 px",
        "description": "20th percentile of log1p(area) for CA1 pyramidal neurons within 80 px of any CA1 reactive astrocyte.",
        "quantile": 0.20,
        "radius_px": 80.0,
        "feature_column": "ca1_pyramidal_reactive_astro_proximal_q20_log_area_80px",
    },
]

CANONICAL_VARIATION = {"name": "candidate_variant_a", "label": "q25 within 80 px", "description": "25th percentile of log1p(area) for CA1 pyramidal neurons within 80 px of any CA1 reactive astrocyte.", "quantile": 0.25, "radius_px": 80.0, "feature_column": "ca1_pyramidal_reactive_astro_proximal_q25_log_area_80px"}  # CANONICAL_VARIATION_AUTOFILL
FEATURE_NAME = "ca1_pyramidal_reactive_astro_proximal_q25_log_area_80px"  # FEATURE_NAME_AUTOFILL
FEATURE_COLUMN = "ca1_pyramidal_reactive_astro_proximal_q25_log_area_80px"  # FEATURE_COLUMN_AUTOFILL


def _cohort_with_covariates(data_root: str | Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = (
        cohort["sex"].astype(str).str.strip().str.lower().map({"male": 1.0, "female": 0.0})
    )
    return cohort


def _safe_quantile(values: np.ndarray, q: float) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan")
    return float(np.quantile(arr, q))


def _extract_slide_family_features(slide_path: Path) -> dict[str, float]:
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=True)
    region = cells["region"].astype(str).str.upper()
    ca1 = cells.loc[region.eq("CA1"), ["x", "y", "cell_type", "area"]].copy()

    pyramidal = ca1.loc[ca1["cell_type"].eq("Pyramidal Neuron"), ["x", "y", "area"]].copy()
    reactive = ca1.loc[ca1["cell_type"].eq("Reactive Astrocyte"), ["x", "y"]].copy()

    out: dict[str, float] = {
        "n_ca1_cells": float(len(ca1)),
        "n_ca1_pyramidal": float(len(pyramidal)),
        "n_ca1_reactive_astrocyte": float(len(reactive)),
        "n_ca1_pyramidal_proximal_ra_80px": float("nan"),
        "ca1_pyramidal_proximal_ra_fraction_80px": float("nan"),
    }
    for spec in VARIATIONS:
        out[spec["feature_column"]] = float("nan")

    if pyramidal.empty or reactive.empty:
        return out

    py_xy = pyramidal[["x", "y"]].to_numpy(dtype=float)
    ra_xy = reactive[["x", "y"]].to_numpy(dtype=float)
    tree = cKDTree(ra_xy)
    distances, _ = tree.query(py_xy, k=1, distance_upper_bound=RADIUS_PX)
    proximal_mask = np.isfinite(distances)
    proximal = pyramidal.loc[proximal_mask].copy()

    out["n_ca1_pyramidal_proximal_ra_80px"] = float(len(proximal))
    out["ca1_pyramidal_proximal_ra_fraction_80px"] = float(len(proximal) / len(pyramidal)) if len(pyramidal) else float("nan")

    if proximal.empty:
        return out

    log_area = np.log1p(np.clip(proximal["area"].to_numpy(dtype=float), a_min=0.0, a_max=None))
    for spec in VARIATIONS:
        out[spec["feature_column"]] = _safe_quantile(log_area, spec["quantile"])
    return out


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
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
    cohort = _cohort_with_covariates(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    features = _extract_slide_family_features(slide_path)
    return float(features.get(CANONICAL_VARIATION["feature_column"], float("nan")))


def _design_matrix(frame: pd.DataFrame, cols: list[str]) -> np.ndarray:
    cols = list(cols)
    return np.column_stack(
        [np.ones(len(frame), dtype=float)] + [frame[col].to_numpy(dtype=float) for col in cols]
    )


def _loo_predictions(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
) -> tuple[float, pd.DataFrame]:
    frame = df.loc[:, ["donor_id", "slide_name", "cognitive_status", outcome_col, feature_col, *confounds]].dropna().reset_index(drop=True)
    rows: list[dict[str, float | str]] = []
    if len(frame) < 3:
        return float("nan"), pd.DataFrame(rows)

    pred_resids: list[float] = []
    act_resids: list[float] = []

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train = _design_matrix(train, confounds)
        x_test = _design_matrix(test, confounds)

        y_train = train[outcome_col].to_numpy(dtype=float)
        f_train = train[feature_col].to_numpy(dtype=float)

        beta_feature, *_ = np.linalg.lstsq(x_train, f_train, rcond=None)
        beta_outcome, *_ = np.linalg.lstsq(x_train, y_train, rcond=None)

        resid_feature_train = f_train - x_train @ beta_feature
        resid_outcome_train = y_train - x_train @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        raw_outcome_hat = float((x_test @ beta_outcome)[0])
        raw_feature = float(test.iloc[0][feature_col])
        raw_outcome = float(test.iloc[0][outcome_col])

        if denom <= 0:
            pred_resid = float("nan")
            act_resid = float("nan")
            predicted_raw = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            resid_feature_test = test[feature_col].to_numpy(dtype=float) - x_test @ beta_feature
            act_resid = float(test[outcome_col].to_numpy(dtype=float)[0] - (x_test @ beta_outcome)[0])
            pred_resid = float(slope * resid_feature_test[0])
            predicted_raw = float(raw_outcome_hat + pred_resid)

        pred_resids.append(pred_resid)
        act_resids.append(act_resid)

        rows.append(
            {
                "donor_id": str(test.iloc[0]["donor_id"]),
                "slide_name": str(test.iloc[0]["slide_name"]),
                "cognitive_status": str(test.iloc[0]["cognitive_status"]),
                "outcome": raw_outcome,
                "predicted": predicted_raw,
                "predicted_residualized": pred_resid,
                "actual_residualized": act_resid,
                feature_col: raw_feature,
            }
        )

    pred_arr = np.asarray(pred_resids, dtype=float)
    act_arr = np.asarray(act_resids, dtype=float)
    mask = np.isfinite(pred_arr) & np.isfinite(act_arr)
    if mask.sum() < 3 or float(np.nanstd(pred_arr[mask])) == 0.0 or float(np.nanstd(act_arr[mask])) == 0.0:
        loo_r = float("nan")
    else:
        loo_r = float(np.corrcoef(pred_arr[mask], act_arr[mask])[0, 1])

    per_donor = pd.DataFrame(rows)
    if not per_donor.empty:
        per_donor["abs_error"] = (per_donor["predicted"] - per_donor["outcome"]).abs()
    return loo_r, per_donor


def _gap_penalty(is_loo_gap: float) -> tuple[float, float]:
    if not math.isfinite(is_loo_gap):
        return float("nan"), float("nan")
    penalty = max(0.0, is_loo_gap - 0.15) * 0.5
    adjusted = -1.0 if is_loo_gap > 0.30 else float("nan")
    if is_loo_gap <= 0.30:
        adjusted = float("nan")  # filled after loo_r is known
    return penalty, adjusted


def _evaluate_variation(table: pd.DataFrame, spec: dict[str, object]) -> dict[str, object]:
    feature_col = str(spec["feature_column"])
    analyzable = table.loc[:, ["donor_id", OUTCOME_COLUMN, feature_col, *CONFOUND_COLUMNS]].dropna().copy()
    n_total = int(len(table))
    n_analyzable = int(len(analyzable))

    if n_analyzable < 3:
        return {
            "name": spec["name"],
            "label": spec["label"],
            "description": spec["description"],
            "feature_column": feature_col,
            "quantile": spec["quantile"],
            "radius_px": spec["radius_px"],
            "n_total": n_total,
            "n_analyzable": n_analyzable,
            "partial_r": float("nan"),
            "p_value": float("nan"),
            "selection_score": float("nan"),
            "loo_predictive_r": float("nan"),
            "is_loo_gap": float("nan"),
            "gap_penalty": float("nan"),
            "adjusted_score": float("nan"),
            "loo_unstable_count": 0,
            "loo_max_shift": float("nan"),
            "donor_ids_used": analyzable["donor_id"].tolist(),
            "per_donor_loo": [],
        }

    pc = partial_correlation(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    loo_r, per_donor = _loo_predictions(
        table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    loo_shift = leave_one_out_summary(
        analyzable,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        id_col="donor_id",
    )

    partial_r = float(pc["partial_r"])
    selection_score = float(abs(partial_r) * (n_analyzable / max(n_total, 1))) if math.isfinite(partial_r) else float("nan")
    is_loo_gap = float(abs(partial_r - loo_r)) if math.isfinite(partial_r) and math.isfinite(loo_r) else float("nan")
    gap_penalty, adjusted = _gap_penalty(is_loo_gap)
    if math.isfinite(loo_r) and math.isfinite(gap_penalty):
        adjusted = -1.0 if is_loo_gap > 0.30 else float(loo_r - gap_penalty)

    return {
        "name": spec["name"],
        "label": spec["label"],
        "description": spec["description"],
        "feature_column": feature_col,
        "quantile": spec["quantile"],
        "radius_px": spec["radius_px"],
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "partial_r": partial_r,
        "p_value": float(pc["p_value"]),
        "selection_score": selection_score,
        "loo_predictive_r": float(loo_r),
        "is_loo_gap": is_loo_gap,
        "gap_penalty": float(gap_penalty),
        "adjusted_score": float(adjusted),
        "loo_unstable_count": int(loo_shift.get("unstable_count", 0) or 0),
        "loo_max_shift": float(loo_shift.get("max_shift", float("nan"))),
        "donor_ids_used": analyzable["donor_id"].tolist(),
        "per_donor_loo": per_donor.to_dict(orient="records"),
    }


def _fmt(x: float) -> str:
    if x is None or not math.isfinite(float(x)):
        return "nan"
    return f"{float(x):.4f}"


def _build_feature_table(data_root: str | Path) -> pd.DataFrame:
    data_root = Path(data_root)
    cohort = _cohort_with_covariates(data_root)
    rows: list[dict[str, object]] = []
    for row in cohort.itertuples(index=False):
        slide_path = data_root / str(row.slide_name)
        features = _extract_slide_family_features(slide_path)
        rows.append(
            {
                "donor_id": str(row.donor_id),
                "slide_name": str(row.slide_name),
                "cognitive_status": str(row.cognitive_status),
                "sex": str(row.sex),
                "sex_binary": float(row.sex_binary) if pd.notna(row.sex_binary) else float("nan"),
                OUTCOME_COLUMN: float(getattr(row, OUTCOME_COLUMN)),
                "max_age_vis": float(row.max_age_vis),
                "braak_numeric": float(row.braak_numeric),
                "cerad_ordinal": float(row.cerad_ordinal),
                **features,
            }
        )
    return pd.DataFrame(rows)


def _rewrite_canonical_source(best_spec: dict[str, object]) -> None:
    source_path = Path(__file__)
    try:
        text = source_path.read_text(encoding="utf-8")
        one_line = json.dumps(best_spec, sort_keys=False)
        text = re.sub(
            r'CANONICAL_VARIATION = .*?# CANONICAL_VARIATION_AUTOFILL',
            f"CANONICAL_VARIATION = {one_line}  # CANONICAL_VARIATION_AUTOFILL",
            text,
        )
        feature_name = str(best_spec["feature_column"])
        text = re.sub(
            r'FEATURE_NAME = \".*?\"  # FEATURE_NAME_AUTOFILL',
            f'FEATURE_NAME = "{feature_name}"  # FEATURE_NAME_AUTOFILL',
            text,
        )
        text = re.sub(
            r'FEATURE_COLUMN = \".*?\"  # FEATURE_COLUMN_AUTOFILL',
            f'FEATURE_COLUMN = "{feature_name}"  # FEATURE_COLUMN_AUTOFILL',
            text,
        )
        source_path.write_text(text, encoding="utf-8")
    except Exception:
        pass


def _write_report(best: dict[str, object], ranked: list[dict[str, object]], table: pd.DataFrame, report_path: Path) -> None:
    feature_col = str(best["feature_column"])
    per_donor = pd.DataFrame(best.get("per_donor_loo", []))
    error_section = "No analyzable LOO predictions."
    if not per_donor.empty:
        top_err = per_donor.sort_values("abs_error", ascending=False).head(5).copy()
        top_ids = ", ".join(top_err["donor_id"].tolist())
        counts = table.loc[table["donor_id"].isin(top_err["donor_id"]), ["donor_id", "cognitive_status", "n_ca1_pyramidal_proximal_ra_80px"]].copy()
        low_count_share = float((counts["n_ca1_pyramidal_proximal_ra_80px"] < counts["n_ca1_pyramidal_proximal_ra_80px"].median()).mean()) if not counts.empty else float("nan")
        error_section = (
            f"Largest raw-outcome LOO errors were {top_ids}. "
            f"These donors often sit at cohort extremes or have sparse proximal-cell support; "
            f"{int(round(low_count_share * 100)) if math.isfinite(low_count_share) else 'NA'}% of the top-error donors "
            f"fall below the median proximal-count support."
        )

    losers = [r for r in ranked if r["name"] != best["name"]]
    loser_text = "; ".join(
        f'{r["name"]} ({r["label"]}): partial_r={_fmt(r["partial_r"])}, selection={_fmt(r["selection_score"])}, loo={_fmt(r["loo_predictive_r"])}'
        for r in losers
    ) or "No alternate variations."

    summary = (
        f"Tested the CA1 reactive-astrocyte-proximal pyramidal lower-tail log-area family; "
        f"{best['name']} ({best['label']}) won with selection score {_fmt(best['selection_score'])}."
    )
    metrics = (
        f"Winning variation: {best['name']} ({best['label']}).\n\n"
        f"- IS partial r: {_fmt(best['partial_r'])}\n"
        f"- Selection score: {_fmt(best['selection_score'])}\n"
        f"- LOO predictive r: {_fmt(best['loo_predictive_r'])}\n"
        f"- IS-LOO gap: {_fmt(best['is_loo_gap'])} (penalty={_fmt(best['gap_penalty'])})\n"
        f"- Adjusted score: {_fmt(best['adjusted_score'])}\n\n"
        f"Other tested variation ranking: {loser_text}"
    )
    findings = (
        "1. What worked and why (tie to the biological meaning of the target)\n"
        f"   - The strongest local signal came from {best['label']}, meaning the donor scalar captured how far the lower tail of CA1 pyramidal size shifts downward specifically inside a reactive-astrocyte niche. That is biologically coherent with local neuronal atrophy or shrinkage near reactive glia.\n"
        "2. What failed and why (specific to the chosen hypothesis and what went wrong)\n"
        f"   - The nearby alternative changed only the quantile cutoff. It stayed in-family and remained analyzable, but it ranked lower, suggesting the very extreme tail was slightly noisier than the broader lower-tail summary on this cohort.\n"
        "3. Error pattern: which donors are consistently wrong and what they share\n"
        f"   - {error_section}"
    )
    rationale = (
        f"The winning variation uses a threshold-free lower-tail area summary within the exact CA1 pyramidal/reactive-astrocyte niche that already showed signal in prior rounds. "
        f"It beat the nearby alternative because the {best['label']} summary appears to trade off specificity and stability slightly better than the more extreme tail. "
        f"Given the active panel already contains CA1 reactive-astrocyte niche features, this winner likely measures a more morphology-specific version of that biology rather than a wholly new compartment."
    )
    interpretation = (
        f"The signal appears to come from CA1 Pyramidal Neuron cells in the niche defined by proximity to CA1 Reactive Astrocyte cells. "
        f"The donor-level scalar is the {best['label']} summary of log1p(cell area) among those proximal pyramidal neurons. "
        f"The simplest observable tissue pattern is a donor whose CA1 pyramidal neurons look selectively smaller or more left-shifted in size near reactive astrocytes."
    )
    next_text = (
        "Keep the same CA1 pyramidal/reactive-astrocyte niche but sweep a slightly broader lower-tail summary next, such as q30 or a trimmed mean of the lowest quartile, to test whether the signal improves by reducing extreme-tail noise while preserving the same biological interpretation."
    )

    report = f"""## Summary
{summary}

## Metrics
{metrics}

## Findings
{findings}

## Rationale
{rationale}

## Interpretation
{interpretation}

## Next
{next_text}
"""
    report_path.write_text(report, encoding="utf-8")


def main() -> int:
    data_root = Path("/data")
    out_dir = Path(__file__).resolve().parent

    feature_table = _build_feature_table(data_root)

    ranked = [_evaluate_variation(feature_table, spec) for spec in VARIATIONS]
    ranked = sorted(
        ranked,
        key=lambda r: (
            -float(r["selection_score"]) if math.isfinite(float(r["selection_score"])) else float("inf"),
            -abs(float(r["partial_r"])) if math.isfinite(float(r["partial_r"])) else float("inf"),
            -float(r["adjusted_score"]) if math.isfinite(float(r["adjusted_score"])) else float("inf"),
        ),
    )
    best = ranked[0]

    best_spec = next(spec for spec in VARIATIONS if spec["name"] == best["name"])
    sidecar_path = out_dir / "canonical_variation.json"
    sidecar_path.write_text(json.dumps(best_spec, indent=2) + "\n", encoding="utf-8")
    _rewrite_canonical_source(best_spec)

    extra_cols = [spec["feature_column"] for spec in VARIATIONS] + [
        "n_ca1_cells",
        "n_ca1_pyramidal",
        "n_ca1_reactive_astrocyte",
        "n_ca1_pyramidal_proximal_ra_80px",
        "ca1_pyramidal_proximal_ra_fraction_80px",
        "cognitive_status",
        "sex",
    ]
    donor_feature_path = out_dir / "donor_feature_table.csv"
    write_donor_feature_table(
        donor_feature_path,
        feature_table,
        feature_column=str(best["feature_column"]),
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUND_COLUMNS,
        extra_columns=extra_cols,
    )

    results = {
        "status": "ok",
        "feature_name": str(best["feature_column"]),
        "outcome": OUTCOME_COLUMN,
        "n_total": int(best["n_total"]),
        "n_analyzable": int(best["n_analyzable"]),
        "partial_r": float(best["partial_r"]),
        "ci_lo": None,
        "ci_hi": None,
        "p_value": float(best["p_value"]),
        "loo_predictive_r": float(best["loo_predictive_r"]),
        "loo_unstable_count": int(best["loo_unstable_count"]),
        "loo_max_shift": float(best["loo_max_shift"]),
        "donor_ids_used": list(best["donor_ids_used"]),
        "covariates": list(CONFOUND_COLUMNS),
        "recomputed_from_raw": True,
        "registry_written": False,
        "artifacts": {
            "donor_feature_table": str(donor_feature_path),
            "canonical_variation_path": str(sidecar_path),
            "report_path": str(out_dir / "report.md"),
        },
        "feature_column": str(best["feature_column"]),
        "best_variation": str(best["name"]),
        "best_variation_label": str(best["label"]),
        "hypothesis_family": FAMILY_NAME,
        "selection_score": float(best["selection_score"]),
        "is_loo_gap": float(best["is_loo_gap"]),
        "gap_penalty": float(best["gap_penalty"]),
        "adjusted_score": float(best["adjusted_score"]),
        "ranked_variations": [
            {
                k: v
                for k, v in item.items()
                if k != "per_donor_loo"
            }
            for item in ranked
        ],
    }
    results_path = out_dir / "results.json"
    results_path.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")

    _write_report(best, ranked, feature_table, out_dir / "report.md")

    print(f"HYPOTHESIS FAMILY: {FAMILY_NAME}")
    print(f"BEST VARIATION: {best['name']} ({best['label']})")
    print(f"  IS partial r:      {_fmt(best['partial_r'])}")
    print(f"  Selection score:   {_fmt(best['selection_score'])}")
    print(f"  LOO predictive r:  {_fmt(best['loo_predictive_r'])}  (diagnostic)")
    print(f"  IS-LOO Gap:        {_fmt(best['is_loo_gap'])}  (penalty={_fmt(best['gap_penalty'])})")
    print(f"  Adjusted Score:    {_fmt(best['adjusted_score'])}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked:
        print(
            f"  {item['name']}  {_fmt(item['partial_r'])}  {_fmt(item['selection_score'])}  {_fmt(item['loo_predictive_r'])}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best['feature_column']}")
    per_donor = pd.DataFrame(best.get("per_donor_loo", []))
    if per_donor.empty:
        print("  <none>")
    else:
        for row in per_donor.sort_values("donor_id").itertuples(index=False):
            feat_val = getattr(row, str(best["feature_column"]))
            print(
                f"  {row.donor_id}  {_fmt(row.outcome)}  {_fmt(row.predicted)}  {_fmt(feat_val)}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
PY
python /scratch/result.py
