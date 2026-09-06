## Summary
Tested the `ca1_reactive_astrocyte_enrichment` family in CA1; `reactive_over_astroglia` won the local sweep with selection score 0.4914.

## Metrics
Winning variation: `reactive_over_astroglia` (`ca1_reactive_astrocyte_enrichment__reactive_over_astroglia`)

- IS partial r: -0.4914
- Selection score: 0.4914
- LOO predictive r: 0.3598
- IS-LOO gap: 0.1316
- Penalty: 0.0000
- Adjusted score: 0.3598
- Bootstrap 95% CI: [-0.7652, -0.1146]
- n analyzable: 35

Other tested variations:
- `reactive_over_astroglia`: partial r -0.4914, selection score 0.4914, LOO r 0.3598
- `reactive_over_all_cells`: partial r -0.4217, selection score 0.4217, LOO r 0.2896

## Findings
1. What worked and why  
   `reactive_over_astroglia` worked best because it isolates CA1 reactive astrocyte state within the astroglial lineage rather than letting the denominator drift with the whole CA1 cellular mixture. That makes it a cleaner donor-level scalar for local gliosis / injury response.

2. What failed and why  
   The nearby denominator swap `reactive_over_all_cells` ranked lower (selection score 0.4217), suggesting that dividing by all CA1 cells diluted the astroglial-state signal with broader CA1 composition.

3. Error pattern: which donors are consistently wrong and what they share  
   The largest absolute LOO prediction errors were in donors H20.33.018, H21.33.023, H21.33.006. Two of the three largest errors (H20.33.018, H21.33.023, H21.33.006) had above-median reactive_over_astroglia values, so the model tends to over-call decline in some donors with strong CA1 gliosis but less severe memory slope than expected.

## Rationale
The best variation is biologically coherent because it measures a shift from homeostatic `Astrocyte` toward `Reactive Astrocyte` specifically inside CA1, the same region already implicated by the accepted panel feature. That framing matches a local injury-response process rather than generic tissue composition. It beat the nearby alternative because the all-cell denominator mixes astroglial activation with large donor-to-donor variation in neurons, oligodendrocytes, and other CA1 cell classes. Because the accepted panel member tracks CA1 pyramidal depletion, this glial-state ratio is biologically distinct and plausibly additive, though panel-level gain must be confirmed by the evaluator.

## Interpretation
The signal appears to reflect CA1 astroglial activation.  
Population: `Reactive Astrocyte` versus `Astrocyte`  
Niche: CA1 annotated tissue compartment  
Feature summary: donor-level CA1 reactive astrocyte enrichment ratio  
Simplest observable pattern: CA1 fields where the astrocyte pool is shifted toward reactive-state labels rather than baseline astrocytes.

## Next
Keep the CA1 reactive-astrocyte lead but test a tighter CA1 niche next, such as reactive-over-astroglia restricted to the CA1 pyramidal layer border or a CA1 reactive-to-pyramidal normalization, since the all-cell denominator lost signal and the remaining errors suggest that diffuse CA1 burden is still too coarse.
