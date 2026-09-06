## Summary
One sentence: tested the ca1_reactive_astrocyte_microclustering family in CA1 astroglia, the winner was reactive_minus_astrocyte_microcluster_fraction_r70px_k3, and its local winner selection score was 0.5125.

## Metrics
Winning variation: reactive_minus_astrocyte_microcluster_fraction_r70px_k3
- IS partial r: -0.5125
- Selection score: 0.5125
- LOO predictive r: 0.4642
- IS-LOO Gap: 0.0484 (penalty=0.0484)
- Adjusted Score: 0.4642

Other tested variations: reactive_minus_astrocyte_microcluster_fraction_r90px_k3 (selection_score=0.4853, partial_r=-0.4853, loo_r=0.4280)

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The winning score worked best when it isolated focal self-clustering of CA1 reactive astrocytes relative to baseline CA1 astrocytes, which is biologically coherent with patchy local gliosis rather than simple bulk abundance.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The broader 90 px radius diluted focality and moved the score toward more generic astroglial packing.
3. Error pattern: which donors are consistently wrong and what they share
   - Largest raw-LOO misses were H20.33.018 (error=0.1641, braak=6, CERAD=3, reactive_frac=0.0956, astro_frac=0.0197); H21.33.040 (error=0.1082, braak=5, CERAD=3, reactive_frac=0.0118, astro_frac=0.0290); H20.33.037 (error=0.0961, braak=5, CERAD=3, reactive_frac=0.0123, astro_frac=0.0160); H21.33.023 (error=0.0959, braak=4, CERAD=0, reactive_frac=0.0045, astro_frac=0.0279); H20.33.043 (error=0.0924, braak=4, CERAD=1, reactive_frac=0.0190, astro_frac=0.1825). These errors tend to come from donors where the reactive-minus-background clustering signal is strong but does not map cleanly onto observed memory-decline severity after confound adjustment.

## Rationale
The best variation is biologically coherent because it compares reactive-astrocyte patching against the background clustering tendency of ordinary CA1 astrocytes, removing some donor-to-donor variation in overall astroglial density. It beat the nearby alternative because the winning neighborhood scale better matched the focal organization of reactive gliosis in CA1. Relative to the current panel, this candidate is most likely to add new information if reactive-astrocyte spatial patchiness carries signal that is not already captured by bulk enrichment or hypertrophy.

## Interpretation
The signal comes from CA1 astroglia. Specifically, the donor scalar is feature__reactive_minus_astrocyte_microcluster_fraction_r70px_k3, the fraction of CA1 reactive astrocytes that sit in same-type microclusters minus the analogous fraction for baseline CA1 astrocytes. The simplest tissue pattern it corresponds to is patchy, locally aggregated islands of reactive astrocytes in CA1 beyond ordinary astroglial packing.

## Next
Keep the CA1 reactive-astrocyte spatial line but next sweep should condition the microclustering score on local pyramidal context or on reactive-astrocyte hypertrophy, because the biggest errors are donors with strong cluster excess whose decline may depend on whether those clusters are neuron-adjacent rather than merely self-clustered.
