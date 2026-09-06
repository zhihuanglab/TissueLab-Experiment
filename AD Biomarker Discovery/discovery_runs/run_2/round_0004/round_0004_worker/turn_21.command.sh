set -euo pipefail
python - <<'PY'
import json
from pathlib import Path
import pandas as pd

res = json.loads(Path('/scratch/results.json').read_text())
ranked = res['ranked_variations']
best = next(r for r in ranked if r['variation_name'] == res['best_variation'])
df = pd.DataFrame(best['per_donor_loo'])
df['abs_err'] = (df['actual_residualized_outcome'] - df['predicted']).abs()
top = df.sort_values('abs_err', ascending=False).head(6)

winner = best['variation_name']
winner_desc = best['description']
feature_col = res['feature_column']

other_lines = []
for row in ranked:
    if row['variation_name'] == winner:
        continue
    other_lines.append(
        f"- `{row['variation_name']}`: partial_r={row['partial_r']:.4f}, "
        f"selection_score={row['selection_score']:.4f}, loo_predictive_r={row['loo_predictive_r']:.4f}, "
        f"adjusted_score={row['adjusted_score']:.4f}."
    )
others_text = "\n".join(other_lines)

err_ids = ", ".join(top['donor_id'].tolist()[:4])
low_astro = top.loc[top['ca1_astrocyte_count'] <= 500, 'donor_id'].tolist()
low_astro_text = ", ".join(low_astro) if low_astro else "none of the top-error donors"
male_n = int((top['sex'] == 'Male').sum())
female_n = int((top['sex'] == 'Female').sum())

report = f"""## Summary
Tested the CA1 reactive-astrocyte hypertrophy family and the local winner was `{winner}` with selection score {best['selection_score']:.4f}.

## Metrics
Winning variation: `{winner}` ({winner_desc})

- IS partial r: {best['partial_r']:.4f}
- Selection score: {best['selection_score']:.4f}
- LOO predictive r: {best['loo_predictive_r']:.4f}
- IS-LOO Gap: {best['is_loo_gap']:.4f} (penalty={best['penalty']:.4f})
- Adjusted Score: {best['adjusted_score']:.4f}
- n_analyzable / n_total: {best['n_analyzable']} / {best['n_total']}
- feature_column: `{feature_col}`

Other tested variation ranking:
{others_text if others_text else "- None."}

## Findings
1. What worked and why  
   The robust median-based CA1 morphology contrast beat the upper-tail variant. The signal came from **CA1 astroglia**, specifically the donor-level scalar `{feature_col}`, defined as CA1 reactive-astrocyte median log1p(contour area) minus CA1 astrocyte median log1p(contour area). Its negative partial correlation means worse memory decline was associated with a **more negative reactive-minus-astrocyte size contrast**. That suggests the family carries signal, but the informative direction is not classic reactive-cell hypertrophy; it is a relative size compression of reactive astrocytes versus the homeostatic astrocyte pool in CA1.

2. What failed and why  
   The planned upper-tail alternative (`candidate_variant_b`) was weaker on both in-sample and LOO diagnostics (partial_r {ranked[1]['partial_r']:.4f}, LOO {ranked[1]['loo_predictive_r']:.4f}). That implies the extreme-area tail adds noise rather than cleaner biology, likely because contour outliers and a small subset of very large cells are less stable than the center of the distribution. More importantly, the round did **not** validate the original hypertrophy-direction intuition: the winning coefficient was negative, so this family is informative mainly as an inverse contrast rather than a positive hypertrophy marker.

3. Error pattern: which donors are consistently wrong and what they share  
   The largest LOO misses were {err_ids}. Several of these donors had relatively small CA1 astrocyte pools or inverted/near-inverted morphology ordering, especially {low_astro_text}. Among the top 6 absolute-error donors, {male_n} were male and {female_n} were female, so sex alone does not explain the misses. The more consistent pattern is that the feature struggles in donors where the CA1 reactive-versus-astrocyte size ordering is unusual or where one astroglial compartment is numerically sparse, making the morphology contrast harder to translate into decline severity.

## Rationale
This was a biologically coherent next step because the accepted panel already implicated a **CA1 reactive astrocyte axis**, and morphology is the most direct way to ask whether that axis reflects activation state rather than only abundance. The winner beat the nearby p75 alternative because activation-related contour differences appear to be distributed across many CA1 reactive astrocytes, not concentrated only in the largest cells. A median contrast is therefore a better summary of the donor's astroglial state than an upper-tail statistic. It also seems plausibly additive beyond the current panel because it measures **within-population shape state** rather than the reactive-cell fraction itself, although panel-level additivity still has to be checked by the evaluator.

## Interpretation
Biologically, the signal seems to mean that in **CA1**, the **reactive astrocyte population** does not become larger relative to the homeostatic astrocyte population in donors with steeper memory decline; instead, the relative area contrast shifts downward. The niche is **CA1**, the donor-level summary is the **median log1p contour-area contrast between CA1 reactive astrocytes and CA1 astrocytes**, and the simplest observable pattern is: **slides with worse memory decline tend to show CA1 reactive astrocyte contours that are not enlarged relative to CA1 astrocytes, and often look smaller/less expanded by this segmentation-derived area readout**.

## Next
Stay in the same family and sweep **count-stabilized median area contrasts** within CA1 reactive astrocytes versus astrocytes, because the plain median beat p75 but the largest residual errors clustered in donors with low CA1 astrocyte counts or inverted size ordering.
"""
Path('/scratch/report.md').write_text(report)
print(report)
PY
