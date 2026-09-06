set -e
python - <<'PY'
import json
from pathlib import Path
import pandas as pd

results = json.loads(Path('/scratch/results.json').read_text())
ranked = results['ranked_variations']
best = next(x for x in ranked if x['variation_name'] == results['best_variation'])
others = [x for x in ranked if x['variation_name'] != results['best_variation']]

cohort = pd.read_csv('/data/training_cohort.csv')
loo = pd.DataFrame(results['best_variation_loo']).sort_values('absolute_error', ascending=False).reset_index(drop=True)
top = loo.head(5).merge(cohort, on=['donor_id', 'slide_name'], how='left')

if len(top):
    error_examples = "; ".join(
        f"{r.donor_id} (|err|={r.absolute_error:.3f}, feature={getattr(r, best['feature_column']):.3f}, Braak={int(r.braak_numeric)}, CERAD={int(r.cerad_ordinal)}, sex={r.sex})"
        for r in top.itertuples(index=False)
    )
    high_feature_misses = int((top[best['feature_column']] > loo[best['feature_column']].median()).sum())
    low_feature_misses = int((top[best['feature_column']] <= loo[best['feature_column']].median()).sum())
    error_pattern = (
        f"The largest errors were mixed: {high_feature_misses} of the top {len(top)} had above-median niche burden "
        f"and {low_feature_misses} had below-median niche burden, so misses came from both unexpectedly severe decline "
        f"at low niche burden and underprediction among some high-burden donors."
    )
else:
    error_examples = "No analyzable donors."
    error_pattern = "No analyzable donors."

other_summary = "; ".join(
    f"{o['variation_name']} ({int(o['radius_px'])} px): partial_r={o['partial_r']:.4f}, selection_score={o['selection_score']:.4f}, loo_r={o['loo_predictive_r']:.4f}"
    for o in others
) if others else "No alternate local variations were tested."

best_col = best['feature_column']
best_radius = int(best['radius_px'])

report = f"""## Summary
Tested the CA1 perineuronal reactive-astrocyte niche family, and candidate_variant_a won with local selection score {best['selection_score']:.4f}.

## Metrics
Winning variation: {results['best_variation']} at {best_radius} px with partial r {best['partial_r']:.4f}, selection score {best['selection_score']:.4f}, LOO predictive r {best['loo_predictive_r']:.4f}, IS-LOO gap {best['is_loo_gap']:.4f}, penalty {best['penalty']:.4f}, and adjusted score {best['adjusted_score']:.4f}. Coverage was {best['n_analyzable']}/{best['n_total']} donors.
Other tested variations: {other_summary}

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The winner was the donor-level fraction of CA1 pyramidal neurons that had at least one CA1 reactive astrocyte within 80 full-resolution pixels. This sharpened the existing CA1 neuron/reactive-astrocyte lead into a spatially specific perineuronal niche, which is biologically plausible for a memory-linked injury response centered on vulnerable CA1 neurons.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The broader 120 px neighborhood was weaker. That suggests the signal is not simply total CA1 reactive-astrocyte burden; expanding the radius likely diluted the intended perineuronal niche by admitting more diffuse astrocytic reactivity that overlaps the already-accepted global CA1 reactive-astrocyte fraction feature.
3. Error pattern: which donors are consistently wrong and what they share
   - Largest LOO errors: {error_examples}
   - {error_pattern}

## Rationale
The best variation is biologically coherent because it measures whether CA1 pyramidal neurons are embedded in a local reactive-astrocyte neighborhood rather than only counting either population globally. It beat the nearby 120 px alternative because the tighter radius better matches a perineuronal scale and preserves local organization. Relative to the current panel, this candidate plausibly adds spatial information beyond overall CA1 pyramidal abundance and reactive-astrocyte lineage fraction, although the positive-but-modest LOO diagnostic and the sign flip versus in-sample partial r argue for caution about instability.

## Interpretation
The signal seems to mean that donors with worse memory decline have a larger share of CA1 pyramidal neurons sitting inside local reactive-astrocyte niches. Population: CA1 pyramidal neurons with CA1 reactive astrocytes as neighbors. Niche: a perineuronal CA1 neighborhood within 80 pixels. Feature summary: donor-level fraction `{best_col}`. Simplest observable pattern: more CA1 neurons appear locally cuffed, flanked, or closely apposed by reactive astrocytes.

## Next
Run the next local sweep around the winning CA1 niche by testing one tighter radius and one requirement for at least two reactive-astrocyte neighbors, to see whether donor-level errors come from overly permissive single-neighbor hits versus a denser perineuronal astrocytic cuff.
"""
Path('/scratch/report.md').write_text(report)
PY
ls -l /scratch/report.md
sed -n '1,220p' /scratch/report.md
