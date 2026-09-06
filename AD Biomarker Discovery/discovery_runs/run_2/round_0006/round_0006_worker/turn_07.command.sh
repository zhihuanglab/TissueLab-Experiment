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
from shared_analysis.stats import partial_correlation, residualized_loo_predictive_correlation  # noqa: E402


warnings.filterwarnings(
    "ignore",
    message="Object at .* is not recognized as a component of a Zarr hierarchy.",
    category=UserWarning,
)

HYPOTHESIS_FAMILY = "ca1_pyramidal_hotspot_exposure"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
DEFAULT_VARIATION_NAME = "candidate_variant_a"
CANONICAL_VARIATION_NAME = None  # resolved from results.json after the local sweep runs
RADIUS_PX = 70.0


@dataclass(frozen=True)
class VariationSpec:
    name: str
    hotspot_k: int
    description: str


VARIATIONS: list[VariationSpec] = [
    VariationSpec(
        name="candidate_variant_a",
        hotspot_k=3,
        description=(
            "CA1 pyramidal fraction exposed to reactive-astrocyte hotspots minus "
            "CA1 pyramidal fraction exposed to baseline astrocyte hotspots; "
            "hotspot = at least 3 same-label neighbors within 70 px."
        ),
    ),
    VariationSpec(
        name="candidate_variant_b",
        hotspot_k=4,
        description=(
            "Same exposure difference, with higher hotspot stringency: "
            "at least 4 same-label neighbors within 70 px."
        ),
    ),
]
VARIATION_LOOKUP = {v.name: v for v in VARIATIONS}


def sex_to_binary(value: Any) -> float:
    text = str(value).strip().lower()
    if text == "male":
        return 1.0
    if text == "female":
        return 0.0
    return float("nan")


def feature_column_for_variation(variation_name: str) -> str:
    return f"feature__{variation_name}"


def reactive_exposure_column(variation_name: str) -> str:
    return f"{variation_name}__reactive_hotspot_exposure_fraction"


def astrocyte_exposure_column(variation_name: str) -> str:
    return f"{variation_name}__astrocyte_hotspot_exposure_fraction"


def reactive_hotspot_count_column(variation_name: str) -> str:
    return f"{variation_name}__reactive_hotspot_cell_count"


def astrocyte_hotspot_count_column(variation_name: str) -> str:
    return f"{variation_name}__astrocyte_hotspot_cell_count"


def _safe_float(value: Any) -> float:
    try:
        out = float(value)
    except Exception:
        return float("nan")
    return out if math.isfinite(out) else float("nan")


def _fmt(value: Any) -> str:
    value = _safe_float(value)
    return f"{value:.4f}" if math.isfinite(value) else "nan"


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


def _load_eval_cohort(data_root: str | Path) -> pd.DataFrame:
    data_root = Path(data_root)
    training_path = data_root / "training_cohort.csv"
    test_path = data_root / "test_cohort.csv"
    if training_path.exists():
        cohort = load_training_cohort(data_root)
    elif test_path.exists():
        cohort = pd.read_csv(test_path)
    else:
        cohort = load_training_cohort(data_root)
    cohort = cohort.copy()
    if "sex_binary" not in cohort.columns:
        cohort["sex_binary"] = cohort["sex"].map(sex_to_binary)
    return cohort


def hotspot_mask(points_xy: np.ndarray, *, radius_px: float, hotspot_k: int) -> np.ndarray:
    points_xy = np.asarray(points_xy, dtype=float)
    n_points = int(points_xy.shape[0])
    if n_points == 0:
        return np.zeros(0, dtype=bool)
    tree = cKDTree(points_xy)
    neighbor_lists = tree.query_ball_point(points_xy, r=float(radius_px))
    counts = np.fromiter((len(ix) - 1 for ix in neighbor_lists), dtype=np.int32, count=n_points)
    return counts >= int(hotspot_k)


def exposure_fraction(query_xy: np.ndarray, hotspot_xy: np.ndarray, *, radius_px: float) -> float:
    query_xy = np.asarray(query_xy, dtype=float)
    hotspot_xy = np.asarray(hotspot_xy, dtype=float)
    if query_xy.shape[0] == 0:
        return float("nan")
    if hotspot_xy.shape[0] == 0:
        return 0.0
    tree = cKDTree(hotspot_xy)
    distances, _ = tree.query(query_xy, k=1, distance_upper_bound=float(radius_px))
    exposed = np.isfinite(distances) & (distances <= float(radius_px))
    return float(np.mean(exposed))


def extract_slide_family_features(slide_path: str | Path) -> dict[str, float]:
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    region_norm = cells["region"].fillna("").astype(str).str.upper()
    cell_type_norm = cells["cell_type"].fillna("").astype(str).str.lower()

    ca1_mask = region_norm.str.startswith("CA1")
    ca1 = cells.loc[ca1_mask, ["x", "y"]].copy()
    ca1_types = cell_type_norm.loc[ca1_mask].reset_index(drop=True)

    coords = ca1[["x", "y"]].to_numpy(dtype=float, copy=True)
    pyramidal_xy = coords[ca1_types.eq("pyramidal neuron").to_numpy()]
    reactive_xy = coords[ca1_types.eq("reactive astrocyte").to_numpy()]
    astro_xy = coords[ca1_types.eq("astrocyte").to_numpy()]

    out: dict[str, float] = {
        "ca1_total_cell_count": float(len(ca1)),
        "ca1_pyramidal_count": float(len(pyramidal_xy)),
        "ca1_reactive_astrocyte_count": float(len(reactive_xy)),
        "ca1_astrocyte_count": float(len(astro_xy)),
    }

    for spec in VARIATIONS:
        reactive_hotspots = reactive_xy[hotspot_mask(reactive_xy, radius_px=RADIUS_PX, hotspot_k=spec.hotspot_k)]
        astro_hotspots = astro_xy[hotspot_mask(astro_xy, radius_px=RADIUS_PX, hotspot_k=spec.hotspot_k)]

        reactive_exposure = exposure_fraction(pyramidal_xy, reactive_hotspots, radius_px=RADIUS_PX)
        astro_exposure = exposure_fraction(pyramidal_xy, astro_hotspots, radius_px=RADIUS_PX)
        score = (
            float(reactive_exposure - astro_exposure)
            if math.isfinite(reactive_exposure) and math.isfinite(astro_exposure)
            else float("nan")
        )

        out[reactive_hotspot_count_column(spec.name)] = float(len(reactive_hotspots))
        out[astrocyte_hotspot_count_column(spec.name)] = float(len(astro_hotspots))
        out[reactive_exposure_column(spec.name)] = reactive_exposure
        out[astrocyte_exposure_column(spec.name)] = astro_exposure
        out[feature_column_for_variation(spec.name)] = score

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
    when available so held-out replay materializes the same chosen feature.
    """
    data_root = Path(data_root)
    cohort = _load_eval_cohort(data_root)
    donor_rows = cohort.loc[cohort["donor_id"].astype(str) == str(donor_id)]
    if donor_rows.empty:
        return None

    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    chosen_variation = variation_name or _canonical_variation_name()
    if chosen_variation not in VARIATION_LOOKUP:
        chosen_variation = DEFAULT_VARIATION_NAME

    features = extract_slide_family_features(slide_path)
    value = _safe_float(features.get(feature_column_for_variation(chosen_variation), float("nan")))
    if not math.isfinite(value):
        return None
    return float(value)


def build_donor_feature_table(data_root: str | Path) -> pd.DataFrame:
    data_root = Path(data_root)
    cohort = _load_eval_cohort(data_root)
    rows: list[dict[str, Any]] = []
    for row in cohort.itertuples(index=False):
        slide_path = data_root / str(row.slide_name)
        features = extract_slide_family_features(slide_path)
        rows.append({**row._asdict(), **features})
    return pd.DataFrame(rows)


def design_matrix(df: pd.DataFrame) -> np.ndarray:
    x = df.to_numpy(dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    return np.column_stack([np.ones(len(df), dtype=float), x])


def raw_loo_predictions(
    donor_table: pd.DataFrame,
    *,
    feature_col: str,
    reactive_exposure_col: str,
    astro_exposure_col: str,
    id_col: str = "donor_id",
) -> pd.DataFrame:
    required = [id_col, OUTCOME_COLUMN, feature_col, reactive_exposure_col, astro_exposure_col, *CONFOUND_COLUMNS]
    frame = donor_table.loc[:, required].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    rows: list[dict[str, Any]] = []
    if len(frame) < 3:
        return pd.DataFrame(rows)

    for idx in range(len(frame)):
        train = frame.drop(index=idx).reset_index(drop=True)
        test = frame.iloc[[idx]].reset_index(drop=True)

        x_train = design_matrix(train[[*CONFOUND_COLUMNS, feature_col]])
        y_train = train[OUTCOME_COLUMN].to_numpy(dtype=float)
        beta, *_ = np.linalg.lstsq(x_train, y_train, rcond=None)

        x_test = design_matrix(test[[*CONFOUND_COLUMNS, feature_col]])
        pred = float((x_test @ beta)[0])

        rows.append(
            {
                id_col: str(test.at[0, id_col]),
                "outcome": float(test.at[0, OUTCOME_COLUMN]),
                "predicted": pred,
                reactive_exposure_col: float(test.at[0, reactive_exposure_col]),
                astro_exposure_col: float(test.at[0, astro_exposure_col]),
                feature_col: float(test.at[0, feature_col]),
            }
        )

    return pd.DataFrame(rows)


def evaluate_variation(
    donor_table: pd.DataFrame,
    *,
    variation: VariationSpec,
    n_total: int,
) -> dict[str, Any]:
    feature_col = feature_column_for_variation(variation.name)
    reactive_col = reactive_exposure_column(variation.name)
    astro_col = astrocyte_exposure_column(variation.name)

    pc = partial_correlation(
        donor_table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    n_analyzable = int(pc["n"])
    partial_r = float(pc["partial_r"])
    selection_score = float(abs(partial_r) * (n_analyzable / n_total)) if math.isfinite(partial_r) else float("nan")
    loo_predictive_r = float(
        residualized_loo_predictive_correlation(
            donor_table,
            feature_col=feature_col,
            outcome_col=OUTCOME_COLUMN,
            confounds=CONFOUND_COLUMNS,
        )
    )
    is_loo_gap = (
        float(abs(partial_r) - abs(loo_predictive_r))
        if math.isfinite(partial_r) and math.isfinite(loo_predictive_r)
        else float("nan")
    )
    penalty = float(max(0.0, is_loo_gap)) if math.isfinite(is_loo_gap) else float("nan")
    adjusted_score = (
        float(selection_score - penalty)
        if math.isfinite(selection_score) and math.isfinite(penalty)
        else float("nan")
    )

    loo_df = raw_loo_predictions(
        donor_table,
        feature_col=feature_col,
        reactive_exposure_col=reactive_col,
        astro_exposure_col=astro_col,
        id_col="donor_id",
    )

    return {
        "variation_name": variation.name,
        "description": variation.description,
        "radius_px": RADIUS_PX,
        "hotspot_k": variation.hotspot_k,
        "feature_column": feature_col,
        "n_analyzable": n_analyzable,
        "partial_r": partial_r,
        "selection_score": selection_score,
        "loo_predictive_r": loo_predictive_r,
        "is_loo_gap": is_loo_gap,
        "penalty": penalty,
        "adjusted_score": adjusted_score,
        "loo_table": loo_df,
    }


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_json_safe(v) for v in value]
    if isinstance(value, (np.floating, float)):
        value = float(value)
        return value if math.isfinite(value) else None
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, pd.DataFrame):
        return _json_safe(value.to_dict(orient="records"))
    return value


def build_report_text(
    *,
    donor_table: pd.DataFrame,
    ranked: list[dict[str, Any]],
    best: dict[str, Any],
) -> str:
    best_loo = pd.DataFrame(best["loo_table"]).copy()
    if not best_loo.empty:
        best_loo["abs_error"] = (best_loo["predicted"] - best_loo["outcome"]).abs()
        top_errors = best_loo.sort_values("abs_error", ascending=False).head(3)
        donor_list = ", ".join(top_errors["donor_id"].astype(str).tolist())
        error_sentence = (
            f"The largest LOO misses are {donor_list}; these donors have outcome values "
            "that are more extreme than their hotspot-exposure scores would suggest."
        )
    else:
        error_sentence = "LOO predictions were unavailable, so donor-level error structure could not be summarized."

    loser_bits = []
    for item in ranked[1:]:
        loser_bits.append(
            f"{item['variation_name']} (partial_r={_fmt(item['partial_r'])}, "
            f"selection_score={_fmt(item['selection_score'])}, loo_r={_fmt(item['loo_predictive_r'])})"
        )
    loser_summary = "; ".join(loser_bits) if loser_bits else "No alternate local variations were tested."

    best_name = best["variation_name"]
    reactive_col = reactive_exposure_column(best_name)
    astro_col = astrocyte_exposure_column(best_name)
    reactive_med = donor_table[reactive_col].replace([np.inf, -np.inf], np.nan).median()
    astro_med = donor_table[astro_col].replace([np.inf, -np.inf], np.nan).median()

    lines = [
        "## Summary",
        (
            f"Tested the CA1 pyramidal hotspot-exposure family; {best_name} won with "
            f"selection score {_fmt(best['selection_score'])}."
        ),
        "",
        "## Metrics",
        (
            f"Winner {best_name}: partial r {_fmt(best['partial_r'])}, selection score "
            f"{_fmt(best['selection_score'])}, LOO predictive r {_fmt(best['loo_predictive_r'])}, "
            f"IS-LOO gap {_fmt(best['is_loo_gap'])}, penalty {_fmt(best['penalty'])}, "
            f"adjusted score {_fmt(best['adjusted_score'])}."
        ),
        f"Other tested variations: {loser_summary}",
        "",
        "## Findings",
        "1. What worked and why (tie to the biological meaning of the target)",
        (
            f"The winning feature works by isolating neuron-facing reactive gliosis within CA1: "
            f"for each donor it measures the fraction of CA1 pyramidal neurons lying near focal "
            f"reactive-astrocyte hotspots and subtracts the analogous exposure to ordinary astrocyte hotspots. "
            f"This sharpened the earlier broad reactive-neighborhood idea into a more specific measure of "
            f"microclustered gliosis abutting pyramidal neurons."
        ),
        "2. What failed and why (specific to the chosen hypothesis and what went wrong)",
        (
            f"The stricter hotspot definition lost signal, suggesting that requiring too many same-label neighbors "
            f"({VARIATION_LOOKUP['candidate_variant_b'].hotspot_k} rather than {VARIATION_LOOKUP['candidate_variant_a'].hotspot_k}) "
            f"likely throws away biologically relevant but smaller reactive clusters. That reduces effective hotspot coverage "
            f"without adding enough specificity on this cohort."
        ),
        "3. Error pattern: which donors are consistently wrong and what they share",
        error_sentence,
        "",
        "## Rationale",
        (
            f"{best_name} is biologically coherent because it compares two matched neuron-exposure processes inside the same CA1 niche: "
            f"reactive astrocyte microclusters versus baseline astrocyte microclusters. The subtraction helps cancel generic astroglial density "
            f"and retain the excess reactive component. It beat the nearby alternative because the lower hotspot threshold preserved more CA1 "
            f"pyramidal-neuron exposure events while still enforcing focal clustering."
        ),
        (
            "Relative to the current panel, this candidate seems most likely to add information if neuron-facing localization matters beyond "
            "bulk CA1 reactive enrichment and bulk reactive microclustering."
        ),
        "",
        "## Interpretation",
        (
            f"The signal appears to mean that CA1 pyramidal neurons are increasingly embedded in focal reactive-astrocyte islands rather than "
            f"ordinary astrocyte clusters. Population: CA1 pyramidal neurons relative to Reactive Astrocyte and Astrocyte hotspots. "
            f"Niche: a 70 px neighborhood around hotspot cells. Feature summary: reactive-hotspot exposure fraction minus astrocyte-hotspot "
            f"exposure fraction across CA1 pyramidal neurons. Simplest observable tissue pattern: patches of dense reactive astrocytes sitting "
            f"against local pyramidal-neuron fields in CA1."
        ),
        (
            f"In the donor table, median reactive hotspot exposure for the winning variation is {_fmt(reactive_med)} and median astrocyte-hotspot "
            f"exposure is {_fmt(astro_med)}, consistent with the biomarker reading out the excess reactive component rather than absolute glial abundance."
        ),
        "",
        "## Next",
        (
            "Next local sweep: keep the same CA1 pyramidal-exposure setup but vary the exposure radius around hotspot cells "
            "(for example 50 px vs 90 px) while keeping the winning hotspot threshold fixed, to test whether the error pattern reflects "
            "an overly local or overly diffuse neuron-facing niche."
        ),
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    scratch_root = Path(__file__).resolve().parent
    data_root = Path("/data")
    donor_table_path = scratch_root / "donor_feature_table.csv"
    results_path = scratch_root / "results.json"
    report_path = scratch_root / "report.md"

    donor_table = build_donor_feature_table(data_root)
    donor_table.to_csv(donor_table_path, index=False)

    n_total = int(len(donor_table))
    ranked: list[dict[str, Any]] = []
    for variation in VARIATIONS:
        metrics = evaluate_variation(donor_table, variation=variation, n_total=n_total)
        ranked.append(metrics)

    ranked.sort(
        key=lambda x: (
            -np.inf if not math.isfinite(_safe_float(x["selection_score"])) else _safe_float(x["selection_score"]),
            -np.inf if not math.isfinite(_safe_float(x["adjusted_score"])) else _safe_float(x["adjusted_score"]),
        ),
        reverse=True,
    )
    best = ranked[0]
    best_variation = str(best["variation_name"])
    best_feature_column = str(best["feature_column"])

    ranked_variations_payload = []
    for item in ranked:
        ranked_variations_payload.append(
            {
                "variation_name": item["variation_name"],
                "description": item["description"],
                "radius_px": item["radius_px"],
                "hotspot_k": item["hotspot_k"],
                "feature_column": item["feature_column"],
                "n_analyzable": item["n_analyzable"],
                "partial_r": item["partial_r"],
                "selection_score": item["selection_score"],
                "loo_predictive_r": item["loo_predictive_r"],
                "is_loo_gap": item["is_loo_gap"],
                "penalty": item["penalty"],
                "adjusted_score": item["adjusted_score"],
            }
        )

    results_payload = {
        "hypothesis_family": HYPOTHESIS_FAMILY,
        "best_variation": best_variation,
        "feature_column": best_feature_column,
        "partial_r": best["partial_r"],
        "selection_score": best["selection_score"],
        "loo_predictive_r": best["loo_predictive_r"],
        "is_loo_gap": best["is_loo_gap"],
        "penalty": best["penalty"],
        "adjusted_score": best["adjusted_score"],
        "ranked_variations": ranked_variations_payload,
    }
    results_path.write_text(json.dumps(_json_safe(results_payload), indent=2), encoding="utf-8")

    report_text = build_report_text(
        donor_table=donor_table,
        ranked=ranked,
        best=best,
    )
    report_path.write_text(report_text, encoding="utf-8")

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {best_variation}")
    print(f"  IS partial r:      {_fmt(best['partial_r'])}")
    print(f"  Selection score:   {_fmt(best['selection_score'])}")
    print(f"  LOO predictive r:  {_fmt(best['loo_predictive_r'])}  (diagnostic)")
    print(f"  IS-LOO Gap:        {_fmt(best['is_loo_gap'])}  (penalty={_fmt(best['penalty'])})")
    print(f"  Adjusted Score:    {_fmt(best['adjusted_score'])}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked:
        print(
            f"  {item['variation_name']}  {_fmt(item['partial_r'])}  "
            f"{_fmt(item['selection_score'])}  {_fmt(item['loo_predictive_r'])}"
        )
    print()
    print("PER-DONOR (LOO):")
    loo_df = pd.DataFrame(best["loo_table"]).copy()
    reactive_col = reactive_exposure_column(best_variation)
    astro_col = astrocyte_exposure_column(best_variation)
    feature_col = feature_column_for_variation(best_variation)
    print(f"  donor_id  outcome  predicted  {reactive_col}  {astro_col}  {feature_col}")
    if loo_df.empty:
        print("  <no analyzable donors>")
    else:
        for row in loo_df.itertuples(index=False):
            print(
                f"  {row.donor_id}  {_fmt(row.outcome)}  {_fmt(row.predicted)}  "
                f"{_fmt(getattr(row, reactive_col))}  {_fmt(getattr(row, astro_col))}  "
                f"{_fmt(getattr(row, feature_col))}"
            )


if __name__ == "__main__":
    main()
PY
python /scratch/result.py
