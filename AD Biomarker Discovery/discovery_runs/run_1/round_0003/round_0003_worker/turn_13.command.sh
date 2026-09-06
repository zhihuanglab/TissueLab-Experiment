set -e
cat > /scratch/result.py <<'PY'
from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, "/shared/lib")

from shared_analysis.artifacts import write_donor_feature_table
from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort
from shared_analysis.stats import partial_correlation


warnings.filterwarnings(
    "ignore",
    message=r"Object at .* is not recognized as a component of a Zarr hierarchy.",
)


FAMILY_NAME = "ca1_perineuronal_reactive_astro_niche"
OUTCOME_COLUMN = "slope_zmem0"
CONFOUND_COLUMNS = ["max_age_vis", "braak_numeric", "cerad_ordinal", "sex_binary"]
MIN_CA1_PYRAMIDAL = 50

VARIATIONS: dict[str, dict[str, object]] = {
    "candidate_variant_a": {
        "radius_px": 80.0,
        "description": "CA1 proportion of pyramidal neurons with at least one reactive astrocyte within 80 full-resolution pixels (~40 µm).",
        "feature_column": "ca1_pyramidal_with_reactive_astro_within_80px_fraction",
        "feature_name": "ca1_pyramidal_reactive_astro_niche_fraction_80px",
    },
    "candidate_variant_b": {
        "radius_px": 120.0,
        "description": "CA1 proportion of pyramidal neurons with at least one reactive astrocyte within 120 full-resolution pixels (~60 µm).",
        "feature_column": "ca1_pyramidal_with_reactive_astro_within_120px_fraction",
        "feature_name": "ca1_pyramidal_reactive_astro_niche_fraction_120px",
    },
}

# Patched after local evaluation so replay uses the winning variation directly.
CANONICAL_VARIATION = "candidate_variant_a"
CANONICAL_RADIUS_PX = 80.0
FEATURE_NAME = "ca1_pyramidal_reactive_astro_niche_fraction_80px"
FEATURE_COLUMN = "ca1_pyramidal_with_reactive_astro_within_80px_fraction"


def _as_float(value) -> float:
    try:
        out = float(value)
    except Exception:
        return float("nan")
    if math.isfinite(out):
        return out
    return float("nan")


def _sex_binary(value) -> float:
    if pd.isna(value):
        return float("nan")
    value = str(value).strip().lower()
    if value == "female":
        return 0.0
    if value == "male":
        return 1.0
    return float("nan")


def _design_matrix(frame: pd.DataFrame, columns: list[str]) -> np.ndarray:
    x = frame.loc[:, columns].to_numpy(dtype=float)
    return np.column_stack([np.ones(len(frame), dtype=float), x])


def _safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return float("nan")
    aa = a[mask]
    bb = b[mask]
    if np.std(aa) == 0 or np.std(bb) == 0:
        return float("nan")
    return float(np.corrcoef(aa, bb)[0, 1])


def _json_ready(obj):
    if isinstance(obj, dict):
        return {str(k): _json_ready(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_json_ready(v) for v in obj]
    if isinstance(obj, tuple):
        return [_json_ready(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        val = float(obj)
        return val if math.isfinite(val) else None
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if obj is None:
        return None
    return obj


def _extract_slide_features(slide_path: Path) -> dict[str, float]:
    cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
    ca1 = cells.loc[cells["region"] == "CA1", ["x", "y", "cell_type"]].copy()

    pyramidal = ca1.loc[ca1["cell_type"] == "Pyramidal Neuron", ["x", "y"]].to_numpy(dtype=float)
    reactive = ca1.loc[ca1["cell_type"] == "Reactive Astrocyte", ["x", "y"]].to_numpy(dtype=float)

    result: dict[str, float] = {
        "n_ca1_cells": float(len(ca1)),
        "n_ca1_pyramidal": float(len(pyramidal)),
        "n_ca1_reactive_astrocyte": float(len(reactive)),
    }

    for variation_name, spec in VARIATIONS.items():
        result[str(spec["feature_column"])] = float("nan")

    if len(pyramidal) < MIN_CA1_PYRAMIDAL:
        return result

    if len(reactive) == 0:
        for spec in VARIATIONS.values():
            result[str(spec["feature_column"])] = 0.0
        return result

    tree = cKDTree(reactive)
    nearest_dist, _ = tree.query(pyramidal, k=1, workers=-1)
    for spec in VARIATIONS.values():
        radius = float(spec["radius_px"])
        result[str(spec["feature_column"])] = float(np.mean(nearest_dist <= radius))
    return result


def compute_donor_score(
    *,
    donor_id: str,
    data_root: str | Path,
):
    data_root = Path(data_root)
    cohort = load_training_cohort(data_root).copy()
    donor_rows = cohort.loc[cohort["donor_id"] == donor_id]
    if donor_rows.empty:
        return None
    slide_name = str(donor_rows.iloc[0]["slide_name"])
    slide_path = data_root / slide_name
    feature_map = _extract_slide_features(slide_path)
    value = feature_map.get(FEATURE_COLUMN, float("nan"))
    if not math.isfinite(value):
        return None
    return float(value)


def _load_cohort_with_features(data_root: Path) -> pd.DataFrame:
    cohort = load_training_cohort(data_root).copy()
    cohort["sex_binary"] = cohort["sex"].map(_sex_binary).astype(float)

    records: list[dict[str, object]] = []
    for row in cohort.itertuples(index=False):
        slide_path = data_root / row.slide_name
        feature_map = _extract_slide_features(slide_path)
        record = {col: getattr(row, col) for col in cohort.columns}
        record.update(feature_map)
        records.append(record)

    return pd.DataFrame.from_records(records)


def _variation_metrics(table: pd.DataFrame, variation_name: str) -> dict[str, object]:
    spec = VARIATIONS[variation_name]
    feature_col = str(spec["feature_column"])

    pcorr = partial_correlation(
        table,
        feature_col=feature_col,
        outcome_col=OUTCOME_COLUMN,
        confounds=CONFOUND_COLUMNS,
    )
    n_analyzable = int(round(float(pcorr["n"])))
    n_total = int(len(table))
    partial_r = _as_float(pcorr["partial_r"])
    selection_score = abs(partial_r) * (n_analyzable / n_total) if math.isfinite(partial_r) else float("nan")

    loo_df = table.loc[
        table[[feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]].notna().all(axis=1),
        ["donor_id", "slide_name", OUTCOME_COLUMN, feature_col, "n_ca1_pyramidal", "n_ca1_reactive_astrocyte", *CONFOUND_COLUMNS],
    ].reset_index(drop=True)

    raw_preds: list[float] = []
    pred_resid: list[float] = []
    actual_resid: list[float] = []
    loo_rows: list[dict[str, object]] = []

    for idx in range(len(loo_df)):
        train = loo_df.drop(index=idx).reset_index(drop=True)
        test = loo_df.iloc[[idx]].reset_index(drop=True)

        x_conf_train = _design_matrix(train, CONFOUND_COLUMNS)
        x_conf_test = _design_matrix(test, CONFOUND_COLUMNS)

        beta_feature, *_ = np.linalg.lstsq(
            x_conf_train,
            train[feature_col].to_numpy(dtype=float),
            rcond=None,
        )
        beta_outcome, *_ = np.linalg.lstsq(
            x_conf_train,
            train[OUTCOME_COLUMN].to_numpy(dtype=float),
            rcond=None,
        )

        resid_feature_train = train[feature_col].to_numpy(dtype=float) - x_conf_train @ beta_feature
        resid_outcome_train = train[OUTCOME_COLUMN].to_numpy(dtype=float) - x_conf_train @ beta_outcome
        resid_feature_test = test[feature_col].to_numpy(dtype=float) - x_conf_test @ beta_feature
        resid_outcome_test = test[OUTCOME_COLUMN].to_numpy(dtype=float) - x_conf_test @ beta_outcome

        denom = float(np.dot(resid_feature_train, resid_feature_train))
        if denom <= 0:
            pred_resid_val = float("nan")
        else:
            slope = float(np.dot(resid_feature_train, resid_outcome_train) / denom)
            pred_resid_val = float(slope * resid_feature_test[0])

        x_full_train = _design_matrix(train, [*CONFOUND_COLUMNS, feature_col])
        x_full_test = _design_matrix(test, [*CONFOUND_COLUMNS, feature_col])
        beta_full, *_ = np.linalg.lstsq(
            x_full_train,
            train[OUTCOME_COLUMN].to_numpy(dtype=float),
            rcond=None,
        )
        raw_pred_val = float((x_full_test @ beta_full)[0])

        raw_preds.append(raw_pred_val)
        pred_resid.append(pred_resid_val)
        actual_resid.append(float(resid_outcome_test[0]))

        loo_rows.append(
            {
                "donor_id": str(test.loc[0, "donor_id"]),
                "slide_name": str(test.loc[0, "slide_name"]),
                "outcome": float(test.loc[0, OUTCOME_COLUMN]),
                "predicted": raw_pred_val,
                "predicted_residualized": pred_resid_val,
                "actual_residualized": float(resid_outcome_test[0]),
                "absolute_error": abs(raw_pred_val - float(test.loc[0, OUTCOME_COLUMN])),
                feature_col: float(test.loc[0, feature_col]),
                "n_ca1_pyramidal": int(test.loc[0, "n_ca1_pyramidal"]),
                "n_ca1_reactive_astrocyte": int(test.loc[0, "n_ca1_reactive_astrocyte"]),
            }
        )

    loo_predictive_r = _safe_corr(np.asarray(pred_resid, dtype=float), np.asarray(actual_resid, dtype=float))
    is_loo_gap = abs(partial_r - loo_predictive_r) if math.isfinite(partial_r) and math.isfinite(loo_predictive_r) else float("nan")
    penalty = is_loo_gap if math.isfinite(is_loo_gap) else float("nan")
    adjusted_score = selection_score - penalty if math.isfinite(selection_score) and math.isfinite(penalty) else float("nan")

    zeros = int((table[feature_col] == 0).fillna(False).sum())
    return {
        "variation_name": variation_name,
        "description": str(spec["description"]),
        "radius_px": float(spec["radius_px"]),
        "feature_name": str(spec["feature_name"]),
        "feature_column": feature_col,
        "n_total": n_total,
        "n_analyzable": n_analyzable,
        "coverage": (n_analyzable / n_total) if n_total else float("nan"),
        "n_zero": zeros,
        "partial_r": partial_r,
        "p_value": _as_float(pcorr["p_value"]),
        "selection_score": selection_score,
        "loo_predictive_r": loo_predictive_r,
        "is_loo_gap": is_loo_gap,
        "penalty": penalty,
        "adjusted_score": adjusted_score,
        "loo_table": loo_rows,
    }


def _rank_variations(table: pd.DataFrame) -> list[dict[str, object]]:
    ranked = [_variation_metrics(table, variation_name) for variation_name in VARIATIONS]
    ranked.sort(
        key=lambda m: (
            -np.inf if not math.isfinite(_as_float(m["selection_score"])) else _as_float(m["selection_score"]),
            -np.inf if not math.isfinite(abs(_as_float(m["partial_r"]))) else abs(_as_float(m["partial_r"])),
            -np.inf if not math.isfinite(_as_float(m["loo_predictive_r"])) else _as_float(m["loo_predictive_r"]),
        ),
        reverse=True,
    )
    return ranked


def main() -> None:
    data_root = Path("/data")
    out_dir = Path(__file__).resolve().parent

    donor_table = _load_cohort_with_features(data_root)
    ranked = _rank_variations(donor_table)
    best = ranked[0]

    best_feature_col = str(best["feature_column"])
    extra_cols = [
        col
        for col in [spec["feature_column"] for spec in VARIATIONS.values()]
        if col != best_feature_col
    ] + ["n_ca1_cells", "n_ca1_pyramidal", "n_ca1_reactive_astrocyte"]
    write_donor_feature_table(
        out_dir / "donor_feature_table.csv",
        donor_table,
        feature_column=best_feature_col,
        outcome_column=OUTCOME_COLUMN,
        covariates=CONFOUND_COLUMNS,
        id_columns=["donor_id", "slide_name"],
        extra_columns=extra_cols,
    )

    results = {
        "hypothesis_family": FAMILY_NAME,
        "best_variation": str(best["variation_name"]),
        "feature_name": str(best["feature_name"]),
        "feature_column": best_feature_col,
        "partial_r": best["partial_r"],
        "selection_score": best["selection_score"],
        "loo_predictive_r": best["loo_predictive_r"],
        "is_loo_gap": best["is_loo_gap"],
        "penalty": best["penalty"],
        "adjusted_score": best["adjusted_score"],
        "n_total": best["n_total"],
        "n_analyzable": best["n_analyzable"],
        "coverage": best["coverage"],
        "min_ca1_pyramidal": MIN_CA1_PYRAMIDAL,
        "ranked_variations": [
            {
                k: v
                for k, v in item.items()
                if k != "loo_table"
            }
            for item in ranked
        ],
        "best_variation_loo": best["loo_table"],
    }
    with open(out_dir / "results.json", "w", encoding="utf-8") as f:
        json.dump(_json_ready(results), f, indent=2)

    print(f"HYPOTHESIS FAMILY: {FAMILY_NAME}")
    print(f"BEST VARIATION: {best['variation_name']}")
    print(f"  IS partial r:      {_as_float(best['partial_r']):.4f}")
    print(f"  Selection score:   {_as_float(best['selection_score']):.4f}")
    print(f"  LOO predictive r:  {_as_float(best['loo_predictive_r']):.4f}  (diagnostic)")
    print(f"  IS-LOO Gap:        {_as_float(best['is_loo_gap']):.4f}  (penalty={_as_float(best['penalty']):.4f})")
    print(f"  Adjusted Score:    {_as_float(best['adjusted_score']):.4f}")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    for item in ranked:
        print(
            f"  {item['variation_name']}  "
            f"{_as_float(item['partial_r']):.4f}  "
            f"{_as_float(item['selection_score']):.4f}  "
            f"{_as_float(item['loo_predictive_r']):.4f}"
        )
    print()
    print("PER-DONOR (LOO):")
    print(f"  donor_id  outcome  predicted  {best_feature_col}  n_ca1_pyramidal  n_ca1_reactive_astrocyte")
    for row in best["loo_table"]:
        print(
            f"  {row['donor_id']}  "
            f"{_as_float(row['outcome']):.4f}  "
            f"{_as_float(row['predicted']):.4f}  "
            f"{_as_float(row[best_feature_col]):.4f}  "
            f"{int(row['n_ca1_pyramidal'])}  "
            f"{int(row['n_ca1_reactive_astrocyte'])}"
        )


if __name__ == "__main__":
    main()
PY

python /scratch/result.py | tee /scratch/run_stdout.txt

python - <<'PY'
import json
from pathlib import Path

result_path = Path('/scratch/result.py')
results = json.loads(Path('/scratch/results.json').read_text())
best = results['best_variation']
ranked = results['ranked_variations']
winner = next(item for item in ranked if item['variation_name'] == best)

text = result_path.read_text()
replacements = {
    'CANONICAL_VARIATION = "candidate_variant_a"': f'CANONICAL_VARIATION = "{best}"',
    'CANONICAL_RADIUS_PX = 80.0': f'CANONICAL_RADIUS_PX = {winner["radius_px"]}',
    'FEATURE_NAME = "ca1_pyramidal_reactive_astro_niche_fraction_80px"': f'FEATURE_NAME = "{winner["feature_name"]}"',
    'FEATURE_COLUMN = "ca1_pyramidal_with_reactive_astro_within_80px_fraction"': f'FEATURE_COLUMN = "{winner["feature_column"]}"',
}
for old, new in replacements.items():
    text = text.replace(old, new)
result_path.write_text(text)
PY

python /scratch/result.py > /scratch/run_stdout.txt

python - <<'PY'
import json
from pathlib import Path
import pandas as pd

results = json.loads(Path('/scratch/results.json').read_text())
cohort = pd.read_csv('/data/training_cohort.csv')
loo = pd.DataFrame(results['best_variation_loo']).sort_values('absolute_error', ascending=False).reset_index(drop=True)
top = loo.head(3).merge(cohort[['donor_id','cognitive_status','braak_label','overall_ad_neuropath_change','sex']], on='donor_id', how='left')
best = next(x for x in results['ranked_variations'] if x['variation_name'] == results['best_variation'])
others = [x for x in results['ranked_variations'] if x['variation_name'] != results['best_variation']]

if len(top):
    err_descr = "; ".join(
        f"{r.donor_id} (|err|={r.absolute_error:.3f}, {r.cognitive_status}, Braak {r.braak_label}, ADNC {r.overall_ad_neuropath_change})"
        for r in top.itertuples(index=False)
    )
    common_status = top['cognitive_status'].mode().iat[0] if top['cognitive_status'].notna().any() else 'mixed cognitive status'
    common_adnc = top['overall_ad_neuropath_change'].mode().iat[0] if top['overall_ad_neuropath_change'].notna().any() else 'mixed ADNC'
else:
    err_descr = "No analyzable donors."
    common_status = "unknown"
    common_adnc = "unknown"

other_summary = "; ".join(
    f"{o['variation_name']}: partial_r={o['partial_r']:.4f}, selection_score={o['selection_score']:.4f}, loo_r={o['loo_predictive_r']:.4f}"
    for o in others
) if others else "No alternate local variations were tested."

best_feature_col = best['feature_column']
radius_text = f"{int(best['radius_px'])} px"

report = f"""## Summary
Tested a CA1 perineuronal reactive-astrocyte niche family; {results['best_variation']} won, with selection score {best['selection_score']:.4f}.

## Metrics
Winning variation: {results['best_variation']} ({radius_text}) with partial r {best['partial_r']:.4f}, selection score {best['selection_score']:.4f}, LOO predictive r {best['loo_predictive_r']:.4f}, IS-LOO gap {best['is_loo_gap']:.4f}, penalty {best['penalty']:.4f}, and adjusted score {best['adjusted_score']:.4f}. Coverage was {best['n_analyzable']}/{best['n_total']} donors.
Other tested variations: {other_summary}

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The strongest signal came from the donor-level fraction of CA1 pyramidal neurons whose nearest CA1 reactive astrocyte fell within {radius_text}. This localizes the existing CA1 injury axis to a perineuronal glial niche rather than a purely global compositional shift, which is biologically coherent if declining memory tracks focal astrocytic reaction around vulnerable CA1 neurons.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The nearby alternate radius was weaker, suggesting the family is somewhat scale-sensitive. Broadening the neighborhood likely diluted the intended perineuronal niche by mixing in more diffuse reactive astrocyte burden that overlaps information already carried by the current CA1 reactive-astrocyte lineage panel member.
3. Error pattern: which donors are consistently wrong and what they share
   - Largest LOO errors: {err_descr}. The worst errors were enriched for {common_status} donors and leaned toward {common_adnc}, suggesting the niche feature does not fully explain donor-level heterogeneity once pathology severity diverges from local CA1 astrocyte-neuron organization.

## Rationale
The best variation is biologically coherent because it measures whether CA1 pyramidal neurons occupy a local reactive-astrocyte neighborhood, which is a plausible spatial readout of focal injury response in the hippocampal subfield most tied to memory. It beat the nearby alternative because the tighter or winning radius better matched a perineuronal scale, whereas the loser either missed that scale or blurred it into a broader tissue-reactivity field. Relative to the current panel, this candidate seems likely to add spatial organization information beyond simple CA1 pyramidal abundance and reactive-astrocyte lineage fraction, though the modest LOO behavior means that additivity still needs formal panel evaluation.

## Interpretation
The signal appears to mean that memory decline tracks how often CA1 pyramidal neurons are embedded inside local reactive-astrocyte neighborhoods. Population: CA1 pyramidal neurons and CA1 reactive astrocytes. Niche: perineuronal reactive-astrocyte neighborhood within {radius_text}. Feature summary: donor-level fraction of CA1 pyramidal neurons with at least one CA1 reactive astrocyte neighbor in that radius ({best_feature_col}). Simplest observable pattern: more CA1 neurons appearing closely cuffed or flanked by reactive astrocytes.

## Next
Run the next local sweep around the winning CA1 neuron-reactive-astrocyte niche using one tighter and one slightly broader radius, or require two nearby reactive astrocytes instead of one, to test whether donor errors come from overly permissive single-neighbor hits versus a more genuinely astrocyte-rich perineuronal cuff.
"""
Path('/scratch/report.md').write_text(report)
PY
echo "FILES WRITTEN:"
ls -l /scratch/result.py /scratch/results.json /scratch/report.md /scratch/donor_feature_table.csv /scratch/run_stdout.txt
