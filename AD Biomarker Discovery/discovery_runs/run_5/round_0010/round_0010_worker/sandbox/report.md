## Summary
One sentence: tested the CA1 peripyramidal reactive astrocyte lymphocyte-contact area burden family, candidate_variant_a won, and its local winner selection score was 0.2107.

## Metrics
Winning variation: candidate_variant_a (ca1_peripyramidal_ra_lymphocyte_area_burden_20um).
- IS partial r: 0.2107
- Selection score: 0.2107
- LOO predictive r: 0.1485
- IS-LOO Gap: 0.0621 (penalty=0.0621)
- Adjusted Score: 0.1485

Ranking:
- candidate_variant_a: partial_r=0.2107, selection_score=0.2107, loo_predictive_r=0.1485
- candidate_variant_b: partial_r=0.2064, selection_score=0.2064, loo_predictive_r=0.1405

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The winning feature isolates CA1 reactive astrocytes that are both peripyramidal and lymphocyte-associated, then upweights them by hypertrophic contour area. That matches the scientific question: the signal is strongest when inflammatory contact and astrocyte enlargement are coupled in the pyramidal niche.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - candidate_variant_b (25 um lymphocyte radius) ranked below the winner with selection score 0.2064, suggesting that broadening the contact definition diluted CA1 peripyramidal specificity.
3. Error pattern: which donors are consistently wrong and what they share
- H20.33.018: outcome -0.304, predicted -0.098, ca1_peripyramidal_ra_lymphocyte_area_burden_20um=0.00
- H21.33.040: outcome -0.033, predicted -0.155, ca1_peripyramidal_ra_lymphocyte_area_burden_20um=0.00
- H21.33.023: outcome -0.155, predicted -0.039, ca1_peripyramidal_ra_lymphocyte_area_burden_20um=0.00

## Rationale
The winning variation is biologically coherent because it focuses on a narrow CA1 niche already implicated by prior accepted features: reactive astrocytes around pyramidal neurons plus local lymphocyte contact. Using summed reactive-astrocyte area over peripyramidal reactive-astrocyte count converts that niche into a severity-weighted burden, which is more specific than a simple presence/count feature. It beat the nearby alternative by preserving the contact scale that best matched the data. It likely carries partially new information beyond the current panel because it fuses hypertrophy with immune contact instead of measuring them separately, though the panel-level additivity still needs formal evaluator review.

## Interpretation
The signal appears to mean a donor-specific burden of enlarged immune-associated reactive astrocytes in the CA1 pyramidal neighborhood. The population is CA1 Reactive Astrocyte cells, the niche is the 30 um peripyramidal zone around CA1 Pyramidal Neurons with nearby CA1 Lymphocytes, the feature summary is summed contact-positive reactive-astrocyte area divided by peripyramidal reactive-astrocyte count, and the simplest observable tissue pattern is clusters of large reactive astrocytes hugging the CA1 pyramidal layer with adjacent lymphocytes.

## Next
Keep the CA1 peripyramidal reactive-astrocyte/lymphocyte burden family, but in the next local sweep change one thing only: retain the winning lymphocyte radius and test whether normalizing by total peripyramidal reactive-astrocyte area rather than count reduces the largest donor-level prediction errors.
