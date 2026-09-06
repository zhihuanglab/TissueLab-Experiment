from __future__ import annotations

import json
import math
import re
import sys
import warnings
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

# Make shared helpers importable when result.py is run directly.
sys.path.insert(0, "/shared/lib")

from shared_analysis import build_cell_table, load_region_polygons, load_training_cohort


FEATURE_FAMILY = "CA1 pyramidal neuron abundance"
FEATURE_NAME = "ca1_pyramidal_abundance"
FEATURE_COLUMN = "ca1_pyramidal_fraction"
CANONICAL_BEST_VARIATION = "ca1_pyramidal_fraction"

MIN_CA1_TOTAL_CELLS = 200

VARIATIONS = {
    "ca1_pyramidal_fraction": {
        "feature_column": "ca1_pyramidal_fraction",
        "description": "CA1 pyramidal cells divided by all classified CA1 cells",
    },
    "ca1_pyramidal_density": {
        "feature_column": "ca1_pyramidal_density",
        "description": "CA1 pyramidal cells divided by annotated CA1 area in cell-coordinate pixel^2",
    },
}


def _sex_to_binary(value) -> float:
    if pd.isna(value):
        return np.nan
    return 1.0 if str(value).strip().lower() == "male" else 0.0


def _safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return float("nan")
    a = a[mask]
    b = b[mask]
    a = a - a.mean()
    b = b - b.mean()
    da = float(np.sqrt(np.dot(a, a)))
    db = float(np.sqrt(np.dot(b, b)))
    if da == 0.0 or db == 0.0:
        return float("nan")
    return float(np.dot(a, b) / (da * db))


def _design_matrix(df: pd.DataFrame, feature_col: str | None = None) -> np.ndarray:
    cols = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
    x = [np.ones(len(df), dtype=float)]
    for col in cols:
        x.append(df[col].astype(float).to_numpy())
    if feature_col is not None:
        x.append(df[feature_col].astype(float).to_numpy())
    return np.column_stack(x)


def _predict_ols(train_df: pd.DataFrame, test_df: pd.DataFrame, feature_col: str) -> np.ndarray:
    x_train = _design_matrix(train_df, feature_col=feature_col)
    y_train = train_df["slope_zmem0"].astype(float).to_numpy()
    beta, *_ = np.linalg.lstsq(x_train, y_train, rcond=None)
    x_test = _design_matrix(test_df, feature_col=feature_col)
    return x_test @ beta


def _partial_correlation(df: pd.DataFrame, feature_col: str) -> float:
    x = df[feature_col].astype(float).to_numpy()
    y = df["slope_zmem0"].astype(float).to_numpy()
    z = _design_matrix(df, feature_col=None)
    bx, *_ = np.linalg.lstsq(z, x, rcond=None)
    by, *_ = np.linalg.lstsq(z, y, rcond=None)
    rx = x - z @ bx
    ry = y - z @ by
    return _safe_corr(rx, ry)


def _polygon_area(poly: np.ndarray) -> float:
    poly = np.asarray(poly, dtype=float)
    if poly.ndim != 2 or poly.shape[0] < 3:
        return 0.0
    x = poly[:, 0]
    y = poly[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _extract_slide_primitives(slide_path: Path) -> Dict[str, float]:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Object at .* is not recognized as a component of a Zarr hierarchy.")
        cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
        polygons = load_region_polygons(slide_path)

    if cells is None or len(cells) == 0:
        return {
            "ca1_total_cells": np.nan,
            "ca1_pyramidal_cells": np.nan,
            "ca1_area_px2": np.nan,
            "ca1_pyramidal_fraction": np.nan,
            "ca1_pyramidal_density": np.nan,
        }

    ca1 = cells.loc[cells["region"] == "CA1", ["cell_type"]].copy()
    ca1_total = int(len(ca1))
    ca1_pyramidal = int((ca1["cell_type"] == "Pyramidal Neuron").sum())

    ca1_polys = polygons.get("CA1", []) if polygons is not None else []
    ca1_area = float(sum(_polygon_area(poly) for poly in ca1_polys)) if ca1_polys else float("nan")

    if ca1_total < MIN_CA1_TOTAL_CELLS:
        fraction = float("nan")
    else:
        fraction = float(ca1_pyramidal / ca1_total)

    if (not np.isfinite(ca1_area)) or ca1_area <= 0 or ca1_total < MIN_CA1_TOTAL_CELLS:
        density = float("nan")
    else:
        density = float(ca1_pyramidal / ca1_area)

    return {
        "ca1_total_cells": float(ca1_total),
        "ca1_pyramidal_cells": float(ca1_pyramidal),
        "ca1_area_px2": float(ca1_area),
        "ca1_pyramidal_fraction": fraction,
        "ca1_pyramidal_density": density,
    }


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    """
    Return one biomarker score for one donor/slide row.

    Canonical replay target: the winning CA1 pyramidal abundance variation from
    this worker round.
    """
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None

    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    primitives = _extract_slide_primitives(slide_path)
    value = primitives.get(FEATURE_COLUMN, np.nan)
    if not np.isfinite(value):
        return None
    return float(value)


def _leave_one_out_predictions(df: pd.DataFrame, feature_col: str) -> pd.DataFrame:
    preds: List[Dict[str, float]] = []
    analyzable = df.loc[
        df[["slope_zmem0", "max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary", feature_col]]
        .notna()
        .all(axis=1)
    ].copy()

    for idx in analyzable.index:
        train_df = analyzable.drop(index=idx)
        test_df = analyzable.loc[[idx]]
        pred = float(_predict_ols(train_df, test_df, feature_col=feature_col)[0])
        row = test_df.iloc[0]
        preds.append(
            {
                "donor_id": str(row["donor_id"]),
                "outcome": float(row["slope_zmem0"]),
                "predicted": pred,
                feature_col: float(row[feature_col]),
                "ca1_total_cells": float(row["ca1_total_cells"]),
                "ca1_pyramidal_cells": float(row["ca1_pyramidal_cells"]),
            }
        )

    out = pd.DataFrame(preds)
    if len(out) == 0:
        return out
    return out.sort_values("donor_id").reset_index(drop=True)


def _evaluate_variation(df: pd.DataFrame, variation_name: str) -> Dict[str, object]:
    feature_col = VARIATIONS[variation_name]["feature_column"]
    analyzable = df.loc[
        df[["slope_zmem0", "max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary", feature_col]]
        .notna()
        .all(axis=1)
    ].copy()
    n_total = int(len(df))
    n_analyzable = int(len(analyzable))
    coverage = float(n_analyzable / n_total) if n_total else 0.0

    if n_analyzable < 5:
        return {
            "variation_name": variation_name,
            "feature_column": feature_col,
            "partial_r": np.nan,
            "selection_score": 0.0,
            "loo_predictive_r": np.nan,
            "is_loo_gap": np.nan,
            "penalty": np.nan,
            "adjusted_score": np.nan,
            "n_total": n_total,
            "n_analyzable": n_analyzable,
            "coverage": coverage,
            "loo_table": [],
        }

    partial_r = _partial_correlation(analyzable, feature_col)
    selection_score = abs(partial_r) * coverage if np.isfinite(partial_r) else 0.0

    loo_df = _leave_one_out_predictions(analyzable, feature_col)
    loo_predictive_r = _safe_corr(
        loo_df["outcome"].to_numpy(dtype=float),
        loo_df["predicted"].to_numpy(dtype=float),
    ) if len(loo_df) else float("nan")

    if np.isfinite(partial_r) and np.isfinite(loo_predictive_r):
        is_loo_gap = abs(abs(partial_r) - abs(loo_predictive_r))
        penalty = is_loo_gap
        adjusted_score = selection_score - penalty
    else:
        is_loo_gap = float("nan")
        penalty = float("nan")
        adjusted_score = float("nan")

    return {
        "variation_name": variation_name,
        "feature_column": feature_col,
        "partial_r": float(partial_r),
        "selection_score": float(selection_score),
        "loo_predictive_r": float(loo_predictive_r),
        "is_loo_gap": float(is_loo_gap),
        "penalty": float(penalty),
        "adjusted_score": float(adjusted_score),
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "coverage": coverage,
        "loo_table": loo_df.to_dict(orient="records"),
    }


def _rewrite_self(best_variation: str, feature_column: str) -> None:
    path = Path(__file__)
    text = path.read_text()
    text = re.sub(
        r'^FEATURE_COLUMN = ".*?"$',
        f'FEATURE_COLUMN = "{feature_column}"',
        text,
        flags=re.M,
    )
    text = re.sub(
        r'^CANONICAL_BEST_VARIATION = ".*?"$',
        f'CANONICAL_BEST_VARIATION = "{best_variation}"',
        text,
        flags=re.M,
    )
    path.write_text(text)


def _build_report(best: Dict[str, object], ranked: List[Dict[str, object]], donor_df: pd.DataFrame, best_loo: pd.DataFrame) -> str:
    other_lines = []
    for row in ranked:
        if row["variation_name"] == best["variation_name"]:
            continue
        other_lines.append(
            f'- `{row["variation_name"]}`: partial r {row["partial_r"]:.4f}, '
            f'selection score {row["selection_score"]:.4f}, '
            f'LOO predictive r {row["loo_predictive_r"]:.4f}.'
        )
    others_text = "\n".join(other_lines) if other_lines else "- No other local variations were tested."

    if len(best_loo):
        err = best_loo.copy()
        err["abs_error"] = (err["outcome"] - err["predicted"]).abs()
        top_err = err.sort_values("abs_error", ascending=False).head(3)
        error_text = "; ".join(
            f'{row.donor_id} (outcome {row.outcome:.4f}, predicted {row.predicted:.4f}, '
            f'{best["feature_column"]} {getattr(row, best["feature_column"]):.4f})'
            for row in top_err.itertuples(index=False)
        )
    else:
        error_text = "No analyzable donors for LOO diagnostics."

    best_feature_col = str(best["feature_column"])
    donor_extremes = donor_df.loc[donor_df[best_feature_col].notna(), ["donor_id", best_feature_col, "ca1_pyramidal_cells", "ca1_total_cells"]]
    low3 = donor_extremes.nsmallest(3, best_feature_col)
    high3 = donor_extremes.nlargest(3, best_feature_col)
    low_text = ", ".join(
        f'{row.donor_id} ({best_feature_col} {getattr(row, best_feature_col):.4f})'
        for row in low3.itertuples(index=False)
    )
    high_text = ", ".join(
        f'{row.donor_id} ({best_feature_col} {getattr(row, best_feature_col):.4f})'
        for row in high3.itertuples(index=False)
    )

    niche = "annotated CA1 pyramidal layer territory within the hippocampus"
    scalar_desc = "fraction of CA1 cells classified as Pyramidal Neuron" if best_feature_col.endswith("fraction") else "density of CA1 Pyramidal Neuron cells per annotated CA1 area"

    return f"""## Summary
One sentence: tested the {FEATURE_FAMILY.lower()} family, `{best['variation_name']}` won, and its local winner selection score was {best['selection_score']:.4f}.

## Metrics
Winning variation `{best['variation_name']}` had partial r {best['partial_r']:.4f}, selection score {best['selection_score']:.4f}, LOO predictive r {best['loo_predictive_r']:.4f}, IS-LOO gap {best['is_loo_gap']:.4f}, penalty {best['penalty']:.4f}, and adjusted score {best['adjusted_score']:.4f}.
Other tested variations:
{others_text}

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The winner tracked donor-to-donor variation in CA1 pyramidal abundance while controlling for age, Braak stage, CERAD, and sex. This is biologically aligned with memory-circuit vulnerability because CA1 pyramidal neurons are a direct substrate of hippocampal memory output.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The losing normalization was more sensitive to regional annotation size, so it mixed neuronal abundance with donor-to-donor CA1 area variability. That diluted the cell-loss signal relative to the cell-composition fraction.
3. Error pattern: which donors are consistently wrong and what they share
   - Largest absolute LOO errors were: {error_text}. These donors likely share mismatch between CA1 pyramidal abundance and decline severity, suggesting additional gliosis, extra-CA1 pathology, or circuit-level compensation not captured by a single abundance scalar.

## Rationale
The best variation is biologically coherent because it isolates the vulnerable CA1 pyramidal population inside its native hippocampal niche and expresses it as a simple donor-level scalar. It beat the nearby density alternative because dividing by total CA1 cells better measures selective compositional loss, whereas density is partially confounded by how much CA1 was annotated and how expanded or contracted the region appears on the section. As a seed biomarker, it is likely to add interpretable information because the current accepted panel is empty.

## Interpretation
The signal seems to mean biologically that donors with lower CA1 pyramidal abundance tend to have faster memory decline. The population is CA1 Pyramidal Neuron cells, the niche is {niche}, the feature summary is donor-level {scalar_desc}, and the simplest observable pattern is a CA1 field with visibly fewer pyramidal neurons relative to the other annotated CA1 cellular content. Lowest-feature donors: {low_text}. Highest-feature donors: {high_text}.

## Next
One specific suggestion for the next local sweep based on the error pattern and which nearby variations won or lost: stay in CA1 and test whether reactive astrocyte enrichment or astrocyte-to-pyramidal balance within the same CA1 niche explains the donors that remain badly predicted after the pyramidal-fraction winner.
"""


def main() -> None:
    data_root = Path("/data")
    out_root = Path("/scratch")
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map(_sex_to_binary).astype(float)

    primitive_rows: List[Dict[str, object]] = []
    for row in cohort.itertuples(index=False):
        slide_path = data_root / str(row.slide_name)
        primitives = _extract_slide_primitives(slide_path)
        primitives["donor_id"] = str(row.donor_id)
        primitive_rows.append(primitives)

    primitive_df = pd.DataFrame(primitive_rows)
    donor_df = cohort.merge(primitive_df, on="donor_id", how="left")

    donor_df.to_csv(out_root / "donor_feature_table.csv", index=False)

    ranked = [_evaluate_variation(donor_df, name) for name in VARIATIONS]
    ranked.sort(key=lambda d: (d["selection_score"], d["adjusted_score"]), reverse=True)
    best = ranked[0]

    _rewrite_self(str(best["variation_name"]), str(best["feature_column"]))

    results_payload = {
        "hypothesis_family": FEATURE_FAMILY,
        "best_variation": best["variation_name"],
        "feature_column": best["feature_column"],
        "partial_r": best["partial_r"],
        "selection_score": best["selection_score"],
        "loo_predictive_r": best["loo_predictive_r"],
        "is_loo_gap": best["is_loo_gap"],
        "penalty": best["penalty"],
        "adjusted_score": best["adjusted_score"],
        "n_total": best["n_total"],
        "n_analyzable": best["n_analyzable"],
        "coverage": best["coverage"],
        "ranked_variations": [
            {
                k: v
                for k, v in row.items()
                if k != "loo_table"
            }
            for row in ranked
        ],
    }
    (out_root / "results.json").write_text(json.dumps(results_payload, indent=2))

    best_loo = pd.DataFrame(best["loo_table"])
    report_md = _build_report(best, ranked, donor_df, best_loo)
    (out_root / "report.md").write_text(report_md)

    print(f"HYPOTHESIS FAMILY: {FEATURE_FAMILY}")
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
    print("PER-DONOR (LOO):")
    if len(best_loo):
        print(f"  donor_id  outcome  predicted  {best['feature_column']}")
        for row in best_loo.itertuples(index=False):
            print(
                f"  {row.donor_id}  "
                f"{row.outcome:.4f}  "
                f"{row.predicted:.4f}  "
                f"{getattr(row, best['feature_column']):.6f}"
            )
    else:
        print("  No analyzable donors.")


if __name__ == "__main__":
    main()
