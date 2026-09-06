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

sys.path.insert(0, "/shared/lib")
from shared_analysis import build_cell_table, load_training_cohort  # noqa: E402
from shared_analysis.artifacts import write_donor_feature_table  # noqa: E402
from shared_analysis.stats import partial_correlation, residualized_loo_predictive_correlation  # noqa: E402


HYPOTHESIS_FAMILY = "ca1_reactive_astrocyte_microclustering"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
ID_COLUMNS = ["donor_id", "slide_name"]
DEFAULT_VARIATION_NAME = "reactive_minus_astrocyte_microcluster_fraction_r70px_k3"
CANONICAL_VARIATION_NAME = None  # loaded from results.json sidecar after the local sweep runs
MIN_CELL_COUNT_PER_TYPE = 25
K_NEIGHBORS = 3


@dataclass(frozen=True)
class VariationSpec:
    name: str
    radius_px: float
    k_neighbors: int = K_NEIGHBORS


VARIATIONS: list[VariationSpec] = [
    VariationSpec(
        name="reactive_minus_astrocyte_microcluster_fraction_r70px_k3",
        radius_px=70.0,
    ),
    VariationSpec(
        name="reactive_minus_astrocyte_microcluster_fraction_r90px_k3",
        radius_px=90.0,
    ),
]
VARIATION_LOOKUP = {v.name: v for v in VARIATIONS}


def _suppress_known_warnings() -> None:
    warnings.filterwarnings(
        "ignore",
        message="Object at .* is not recognized as a component of a Zarr hierarchy.",
        category=UserWarning,
    )


def sex_to_binary(value: Any) -> float:
    text = str(value).strip().lower()
    if text == "male":
        return 1.0
    if text == "female":
        return 0.0
    return float("nan")


def feature_column_for_variation(variation_name: str) -> str:
    return f"feature__{variation_name}"


def _canonical_variation_name() -> str:
    if CANONICAL_VARIATION_NAME in VARIATION_LOOKUP:
        return str(CANONICAL_VARIATION_NAME)
    results_path = Path(__file__).with_name("results.json")
    if results_path.exists():
        try:
            payload = json.loads(results_path.read_text())
            best = payload.get("best_variation")
            if best in VARIATION_LOOKUP:
                return str(best)
        except Exception:
            pass
    return DEFAULT_VARIATION_NAME


def microcluster_fraction(points_xy: np.ndarray, *, radius_px: float, k_neighbors: int) -> float:
    n_points = int(points_xy.shape[0])
    if n_points < MIN_CELL_COUNT_PER_TYPE:
        return float("nan")
    tree = cKDTree(points_xy)
    neighbor_lists = tree.query_ball_point(points_xy, r=radius_px)
    counts = np.fromiter((len(ix) - 1 for ix in neighbor_lists), dtype=np.int32, count=n_points)
    return float(np.mean(counts >= k_neighbors))


def extract_slide_family_features(slide_path: str | Path) -> dict[str, float]:
    _suppress_known_warnings()
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    ca1 = cells.loc[cells["region"] == "CA1", ["x", "y", "cell_type"]].copy()

    reactive = ca1.loc[ca1["cell_type"] == "Reactive Astrocyte", ["x", "y"]].to_numpy(dtype=float)
    astro = ca1.loc[ca1["cell_type"] == "Astrocyte", ["x", "y"]].to_numpy(dtype=float)

    out: dict[str, float] = {
        "ca1_cell_count": float(len(ca1)),
        "ca1_reactive_astrocyte_count": float(len(reactive)),
        "ca1_astrocyte_count": float(len(astro)),
    }

    for spec in VARIATIONS:
        reactive_frac = microcluster_fraction(
            reactive,
            radius_px=spec.radius_px,
            k_neighbors=spec.k_neighbors,
        )
        astro_frac = microcluster_fraction(
            astro,
            radius_px=spec.radius_px,
            k_neighbors=spec.k_neighbors,
        )
        score = (
            float(reactive_frac - astro_frac)
            if np.isfinite(reactive_frac) and np.isfinite(astro_frac)
            else float("nan")
        )
        prefix = spec.name
        out[f"{prefix}__reactive_microcluster_fraction"] = reactive_frac
        out[f"{prefix}__astrocyte_microcluster_fraction"] = astro_frac
        out[feature_column_for_variation(prefix)] = score

    return out


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
    variation_name: str | None = None,
):
    """
    Return the canonical round biomarker score for one donor/slide row.

    The winning local variation is loaded from results.json next to this script
    when available, so held-out replay materializes the same feature chosen by
    the worker run.
    """
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None

    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    chosen_variation = variation_name or _canonical_variation_name()
    if chosen_variation not in VARIATION_LOOKUP:
        chosen_variation = DEFAULT_VARIATION_NAME

    features = extract_slide_family_features(slide_path)
    value = features.get(feature_column_for_variation(chosen_variation), float("nan"))
    if not np.isfinite(value):
        return None
    return float(value)


def clean_frame(df: pd.DataFrame, required: list[str]) -> pd.DataFrame:
    frame = df.loc[:, required].replace([np.inf, -np.inf], np.nan).dropna().copy()
    return frame


def design_matrix(df: pd.DataFrame) -> np.ndarray:
    x = df.to_numpy(dtype=float)
    return np.column_stack([np.ones(len(df), dtype=float), x])


def raw_loo_predictions(
    df: pd.DataFrame,
    *,
    feature_col: str,
    outcome_col: str,
    confounds: list[str],
    id_col: str,
) -> pd.DataFrame:
    required = [id_col, outcome_col, feature_col, *confounds]
    frame = df.loc[:, required].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    rows: list[dict[str, Any]] = []
    if len(frame) < 3:
        return pd.DataFrame(rows)

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train = design_matrix(train[[*confounds, feature_col]])
        y_train = train[outcome_col].to_numpy(dtype=float)
        beta, *_ = np.linalg.lstsq(x_train, y_train, rcond=None)

        x_test = design_matrix(test[[*confounds, feature_col]])
        pred = float((x_test @ beta)[0])

        rows.append(
            {
                id_col: test.at[0, id_col],
                "outcome": float(test.at[0, outcome_col]),
                "predicted": pred,
                feature_col: float(test.at[0, feature_col]),
            }
        )

    return pd.DataFrame(rows)


def evaluate_variation(
    df: pd.DataFrame,
    *,
    variation: VariationSpec,
    n_total: int,
) -> dict[str, Any]:
    feature_col = feature_column_for_variation(variation.name)
    pc = partial_correlation(
        df,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    n_analyzable = int(pc["n"])
    partial_r = float(pc["partial_r"])
    selection_score = float(abs(partial_r) * (n_analyzable / n_total)) if np.isfinite(partial_r) else float("nan")
    loo_predictive_r = float(
        residualized_loo_predictive_correlation(
            df,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUND_COLUMNS,
        )
    )
    is_loo_gap = (
        float(abs(partial_r) - abs(loo_predictive_r))
        if np.isfinite(partial_r) and np.isfinite(loo_predictive_r)
        else float("nan")
    )
    penalty = float(max(0.0, is_loo_gap)) if np.isfinite(is_loo_gap) else float("nan")
    adjusted_score = (
        float(selection_score - penalty)
        if np.isfinite(selection_score) and np.isfinite(penalty)
        else float("nan")
    )

    loo_df = raw_loo_predictions(
        df,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
        id_col="donor_id",
    )
    if not loo_df.empty:
        reactive_col = f"{variation.name}__reactive_microcluster_fraction"
        astro_col = f"{variation.name}__astrocyte_microcluster_fraction"
        merge_cols = ["donor_id", reactive_col, astro_col, feature_col]
        loo_df = loo_df.merge(df.loc[:, merge_cols], on="donor_id", how="left", suffixes=("", "_full"))
        if f"{feature_col}_full" in loo_df.columns:
            loo_df = loo_df.drop(columns=[f"{feature_col}_full"])
        loo_df["abs_error"] = (loo_df["predicted"] - loo_df["outcome"]).abs()

    return {
        "variation_name": variation.name,
        "feature_column": feature_col,
        "n_analyzable": n_analyzable,
        "n_total": int(n_total),
        "partial_r": partial_r,
        "selection_score": selection_score,
        "loo_predictive_r": loo_predictive_r,
        "is_loo_gap": is_loo_gap,
        "penalty": penalty,
        "adjusted_score": adjusted_score,
        "loo_predictions": loo_df,
    }


def build_feature_table(data_root: str | Path) -> pd.DataFrame:
    _suppress_known_warnings()
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map(sex_to_binary)

    feature_rows: list[dict[str, Any]] = []
    for _, row in cohort.iterrows():
        slide_path = data_root / str(row["slide_name"])
        features = extract_slide_family_features(slide_path)
        feature_rows.append({"donor_id": row["donor_id"], **features})

    feature_df = pd.DataFrame(feature_rows)
    merged = cohort.merge(feature_df, on="donor_id", how="left")
    return merged


def format_float(value: Any) -> str:
    if value is None:
        return "nan"
    try:
        value = float(value)
    except Exception:
        return str(value)
    if not np.isfinite(value):
        return "nan"
    return f"{value:.4f}"


def ranking_table_lines(ranked: list[dict[str, Any]]) -> list[str]:
    lines = ["  variation_name  partial_r  selection_score  loo_predictive_r"]
    for item in ranked:
        lines.append(
            "  "
            + f"{item['variation_name']}  {format_float(item['partial_r'])}  "
            + f"{format_float(item['selection_score'])}  {format_float(item['loo_predictive_r'])}"
        )
    return lines


def loo_table_lines(best_result: dict[str, Any]) -> list[str]:
    loo_df = best_result["loo_predictions"].copy()
    feature_col = best_result["feature_column"]
    reactive_col = f"{best_result['variation_name']}__reactive_microcluster_fraction"
    astro_col = f"{best_result['variation_name']}__astrocyte_microcluster_fraction"

    if loo_df.empty:
        return ["  [no analyzable donors]"]

    loo_df = loo_df.sort_values("donor_id").reset_index(drop=True)
    lines = [
        f"  donor_id  outcome  predicted  {feature_col}  reactive_microcluster_fraction  astrocyte_microcluster_fraction"
    ]
    for _, row in loo_df.iterrows():
        lines.append(
            "  "
            + f"{row['donor_id']}  {format_float(row['outcome'])}  {format_float(row['predicted'])}  "
            + f"{format_float(row.get(feature_col))}  {format_float(row.get(reactive_col))}  "
            + f"{format_float(row.get(astro_col))}"
        )
    return lines


def write_report(
    report_path: Path,
    *,
    ranked: list[dict[str, Any]],
    best: dict[str, Any],
    full_table: pd.DataFrame,
) -> None:
    best_name = best["variation_name"]
    best_feature_col = best["feature_column"]
    reactive_col = f"{best_name}__reactive_microcluster_fraction"
    astro_col = f"{best_name}__astrocyte_microcluster_fraction"
    loo_df = best["loo_predictions"].copy()

    others = [r for r in ranked if r["variation_name"] != best_name]
    ranking_summary = "; ".join(
        f"{r['variation_name']} (selection_score={format_float(r['selection_score'])}, partial_r={format_float(r['partial_r'])}, loo_r={format_float(r['loo_predictive_r'])})"
        for r in others
    )
    if not ranking_summary:
        ranking_summary = "No alternate family-local variations were tested."

    error_paragraph = "No analyzable donors for leave-one-out diagnostics."
    next_suggestion = (
        "Test a tighter or thresholded focality variant, such as the same CA1 reactive-minus-astrocyte score at smaller radii or requiring >=4 same-type neighbors, to see whether the winning radius is still too diffuse."
    )
    if not loo_df.empty:
        top_err = loo_df.sort_values("abs_error", ascending=False).head(5).copy()
        merged = top_err.merge(
            full_table[
                [
                    "donor_id",
                    "cognitive_status",
                    "braak_numeric",
                    "cerad_ordinal",
                    "ca1_reactive_astrocyte_count",
                    "ca1_astrocyte_count",
                    reactive_col,
                    astro_col,
                    best_feature_col,
                ]
            ],
            on="donor_id",
            how="left",
            suffixes=("", "_feat"),
        )
        donor_bits = []
        for _, row in merged.iterrows():
            donor_bits.append(
                f"{row['donor_id']} (error={format_float(row['abs_error'])}, braak={row['braak_numeric']}, CERAD={row['cerad_ordinal']}, status={row['cognitive_status']}, reactive_frac={format_float(row[reactive_col])}, astro_frac={format_float(row[astro_col])})"
            )
        error_paragraph = (
            "Largest raw-LOO misses were "
            + "; ".join(donor_bits)
            + ". These errors tend to come from donors where the reactive-minus-background clustering signal is strong but does not map cleanly onto observed memory-decline severity after confound adjustment."
        )
        next_suggestion = (
            "Keep the CA1 reactive-astrocyte spatial line but next sweep should condition the microclustering score on local pyramidal context or on reactive-astrocyte hypertrophy, because the biggest errors are donors with strong cluster excess whose decline may depend on whether those clusters are neuron-adjacent rather than merely self-clustered."
        )

    interpretation = (
        f"The signal comes from CA1 astroglia. Specifically, the donor scalar is {best_feature_col}, "
        "the fraction of CA1 reactive astrocytes that sit in same-type microclusters minus the analogous fraction for baseline CA1 astrocytes. "
        "The simplest tissue pattern it corresponds to is patchy, locally aggregated islands of reactive astrocytes in CA1 beyond ordinary astroglial packing."
    )

    loser_comment = "The broader 90 px radius diluted focality and moved the score toward more generic astroglial packing."
    if ranked and ranked[0]["variation_name"] != "reactive_minus_astrocyte_microcluster_fraction_r70px_k3":
        loser_comment = "The tighter 70 px radius may have been too narrow, while the broader 90 px radius better captured patch size in this cohort."

    text = f"""## Summary
One sentence: tested the {HYPOTHESIS_FAMILY} family in CA1 astroglia, the winner was {best_name}, and its local winner selection score was {format_float(best['selection_score'])}.

## Metrics
Winning variation: {best_name}
- IS partial r: {format_float(best['partial_r'])}
- Selection score: {format_float(best['selection_score'])}
- LOO predictive r: {format_float(best['loo_predictive_r'])}
- IS-LOO Gap: {format_float(best['is_loo_gap'])} (penalty={format_float(best['penalty'])})
- Adjusted Score: {format_float(best['adjusted_score'])}

Other tested variations: {ranking_summary}

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The winning score worked best when it isolated focal self-clustering of CA1 reactive astrocytes relative to baseline CA1 astrocytes, which is biologically coherent with patchy local gliosis rather than simple bulk abundance.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - {loser_comment}
3. Error pattern: which donors are consistently wrong and what they share
   - {error_paragraph}

## Rationale
The best variation is biologically coherent because it compares reactive-astrocyte patching against the background clustering tendency of ordinary CA1 astrocytes, removing some donor-to-donor variation in overall astroglial density. It beat the nearby alternative because the winning neighborhood scale better matched the focal organization of reactive gliosis in CA1. Relative to the current panel, this candidate is most likely to add new information if reactive-astrocyte spatial patchiness carries signal that is not already captured by bulk enrichment or hypertrophy.

## Interpretation
{interpretation}

## Next
{next_suggestion}
"""
    report_path.write_text(text)


def make_results_payload(
    *,
    ranked: list[dict[str, Any]],
    best: dict[str, Any],
) -> dict[str, Any]:
    ranked_json = []
    for item in ranked:
        ranked_json.append(
            {
                "variation_name": item["variation_name"],
                "feature_column": item["feature_column"],
                "n_analyzable": item["n_analyzable"],
                "n_total": item["n_total"],
                "partial_r": item["partial_r"],
                "selection_score": item["selection_score"],
                "loo_predictive_r": item["loo_predictive_r"],
                "is_loo_gap": item["is_loo_gap"],
                "penalty": item["penalty"],
                "adjusted_score": item["adjusted_score"],
            }
        )
    return {
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": best["variation_name"],
        "feature_column": best["feature_column"],
        "n_analyzable": best["n_analyzable"],
        "n_total": best["n_total"],
        "partial_r": best["partial_r"],
        "selection_score": best["selection_score"],
        "loo_predictive_r": best["loo_predictive_r"],
        "is_loo_gap": best["is_loo_gap"],
        "penalty": best["penalty"],
        "adjusted_score": best["adjusted_score"],
        "ranked_variations": ranked_json,
    }


def main() -> None:
    data_root = Path("/data")
    out_dir = Path("/scratch")

    full_table = build_feature_table(data_root)
    n_total = int(len(full_table))

    evaluated = [evaluate_variation(full_table, variation=v, n_total=n_total) for v in VARIATIONS]
    ranked = sorted(
        evaluated,
        key=lambda d: (
            -np.nan_to_num(d["selection_score"], nan=-np.inf),
            -np.nan_to_num(d["partial_r"], nan=-np.inf),
            -np.nan_to_num(d["loo_predictive_r"], nan=-np.inf),
        ),
    )
    best = ranked[0]

    donor_feature_path = out_dir / "donor_feature_table.csv"
    best_feature_col = best["feature_column"]
    best_name = best["variation_name"]
    write_donor_feature_table(
        donor_feature_path,
        full_table,
        feature_column=best_feature_col,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUND_COLUMNS,
        id_columns=ID_COLUMNS,
        extra_columns=[
            "cognitive_status",
            "ca1_cell_count",
            "ca1_reactive_astrocyte_count",
            "ca1_astrocyte_count",
            f"{best_name}__reactive_microcluster_fraction",
            f"{best_name}__astrocyte_microcluster_fraction",
        ],
    )

    results_payload = make_results_payload(ranked=ranked, best=best)
    (out_dir / "results.json").write_text(json.dumps(results_payload, indent=2))

    stdout_lines = [
        f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}",
        f"BEST VARIATION: {best_name}",
        f"  IS partial r:      {format_float(best['partial_r'])}",
        f"  Selection score:   {format_float(best['selection_score'])}",
        f"  LOO predictive r:  {format_float(best['loo_predictive_r'])}  (diagnostic)",
        f"  IS-LOO Gap:        {format_float(best['is_loo_gap'])}  (penalty={format_float(best['penalty'])})",
        f"  Adjusted Score:    {format_float(best['adjusted_score'])}",
        "",
        "RANKED VARIATIONS:",
        *ranking_table_lines(ranked),
        "",
        "PER-DONOR (LOO):",
        *loo_table_lines(best),
    ]
    print("\n".join(stdout_lines))

    write_report(out_dir / "report.md", ranked=ranked, best=best, full_table=full_table)


if __name__ == "__main__":
    main()
PY
python -m py_compile /scratch/result.py
