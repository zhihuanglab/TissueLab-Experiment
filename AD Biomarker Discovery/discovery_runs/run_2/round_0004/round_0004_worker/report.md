## Summary
Tested the CA1 reactive-astrocyte hypertrophy family and the local winner was `candidate_variant_a` with selection score 0.3576.

## Metrics
Winning variation: `candidate_variant_a` (Baseline CA1 reactive astrocyte hypertrophy contrast: median(log1p(area_px2)) of CA1 Reactive Astrocytes minus median(log1p(area_px2)) of CA1 Astrocytes.)

- IS partial r: -0.3576
- Selection score: 0.3576
- LOO predictive r: 0.2201
- IS-LOO Gap: 0.1375 (penalty=0.0000)
- Adjusted Score: 0.2201
- n_analyzable / n_total: 35 / 35
- feature_column: `reactive_minus_astrocyte_median_log1p_area_ca1`

Other tested variation ranking:
- `candidate_variant_b`: partial_r=-0.3135, selection_score=0.3135, loo_predictive_r=0.1603, adjusted_score=0.1587.

## Findings
1. What worked and why  
   The robust median-based CA1 morphology contrast beat the upper-tail variant. The signal came from **CA1 astroglia**, specifically the donor-level scalar `reactive_minus_astrocyte_median_log1p_area_ca1`, defined as CA1 reactive-astrocyte median log1p(contour area) minus CA1 astrocyte median log1p(contour area). Its negative partial correlation means worse memory decline was associated with a **more negative reactive-minus-astrocyte size contrast**. That suggests the family carries signal, but the informative direction is not classic reactive-cell hypertrophy; it is a relative size compression of reactive astrocytes versus the homeostatic astrocyte pool in CA1.

2. What failed and why  
   The planned upper-tail alternative (`candidate_variant_b`) was weaker on both in-sample and LOO diagnostics (partial_r -0.3135, LOO 0.1603). That implies the extreme-area tail adds noise rather than cleaner biology, likely because contour outliers and a small subset of very large cells are less stable than the center of the distribution. More importantly, the round did **not** validate the original hypertrophy-direction intuition: the winning coefficient was negative, so this family is informative mainly as an inverse contrast rather than a positive hypertrophy marker.

3. Error pattern: which donors are consistently wrong and what they share  
   The largest LOO misses were H20.33.018, H21.33.023, H20.33.038, H21.33.004. Several of these donors had relatively small CA1 astrocyte pools or inverted/near-inverted morphology ordering, especially H20.33.018, H20.33.038. Among the top 6 absolute-error donors, 4 were male and 2 were female, so sex alone does not explain the misses. The more consistent pattern is that the feature struggles in donors where the CA1 reactive-versus-astrocyte size ordering is unusual or where one astroglial compartment is numerically sparse, making the morphology contrast harder to translate into decline severity.

## Rationale
This was a biologically coherent next step because the accepted panel already implicated a **CA1 reactive astrocyte axis**, and morphology is the most direct way to ask whether that axis reflects activation state rather than only abundance. The winner beat the nearby p75 alternative because activation-related contour differences appear to be distributed across many CA1 reactive astrocytes, not concentrated only in the largest cells. A median contrast is therefore a better summary of the donor's astroglial state than an upper-tail statistic. It also seems plausibly additive beyond the current panel because it measures **within-population shape state** rather than the reactive-cell fraction itself, although panel-level additivity still has to be checked by the evaluator.

## Interpretation
Biologically, the signal seems to mean that in **CA1**, the **reactive astrocyte population** does not become larger relative to the homeostatic astrocyte population in donors with steeper memory decline; instead, the relative area contrast shifts downward. The niche is **CA1**, the donor-level summary is the **median log1p contour-area contrast between CA1 reactive astrocytes and CA1 astrocytes**, and the simplest observable pattern is: **slides with worse memory decline tend to show CA1 reactive astrocyte contours that are not enlarged relative to CA1 astrocytes, and often look smaller/less expanded by this segmentation-derived area readout**.

## Next
Stay in the same family and sweep **count-stabilized median area contrasts** within CA1 reactive astrocytes versus astrocytes, because the plain median beat p75 but the largest residual errors clustered in donors with low CA1 astrocyte counts or inverted size ordering.
