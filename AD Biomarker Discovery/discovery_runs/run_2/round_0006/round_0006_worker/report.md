## Summary
Tested the CA1 pyramidal hotspot-exposure family; candidate_variant_a won with selection score 0.5927.

## Metrics
Winner candidate_variant_a: partial r -0.5927, selection score 0.5927, LOO predictive r 0.5431, IS-LOO gap 0.0496, penalty 0.0496, adjusted score 0.5431.
Other tested variations: candidate_variant_b (partial_r=-0.5269, selection_score=0.5269, loo_r=0.4464)

## Findings
1. What worked and why (tie to the biological meaning of the target)
The winning feature works by isolating neuron-facing reactive gliosis within CA1: for each donor it measures the fraction of CA1 pyramidal neurons lying near focal reactive-astrocyte hotspots and subtracts the analogous exposure to ordinary astrocyte hotspots. This sharpened the earlier broad reactive-neighborhood idea into a more specific measure of microclustered gliosis abutting pyramidal neurons.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
The stricter hotspot definition lost signal, suggesting that requiring too many same-label neighbors (4 rather than 3) likely throws away biologically relevant but smaller reactive clusters. That reduces effective hotspot coverage without adding enough specificity on this cohort.
3. Error pattern: which donors are consistently wrong and what they share
The largest LOO misses are H20.33.018, H20.33.046, H20.33.037; these donors have outcome values that are more extreme than their hotspot-exposure scores would suggest.

## Rationale
candidate_variant_a is biologically coherent because it compares two matched neuron-exposure processes inside the same CA1 niche: reactive astrocyte microclusters versus baseline astrocyte microclusters. The subtraction helps cancel generic astroglial density and retain the excess reactive component. It beat the nearby alternative because the lower hotspot threshold preserved more CA1 pyramidal-neuron exposure events while still enforcing focal clustering.
Relative to the current panel, this candidate seems most likely to add information if neuron-facing localization matters beyond bulk CA1 reactive enrichment and bulk reactive microclustering.

## Interpretation
The signal appears to mean that CA1 pyramidal neurons are increasingly embedded in focal reactive-astrocyte islands rather than ordinary astrocyte clusters. Population: CA1 pyramidal neurons relative to Reactive Astrocyte and Astrocyte hotspots. Niche: a 70 px neighborhood around hotspot cells. Feature summary: reactive-hotspot exposure fraction minus astrocyte-hotspot exposure fraction across CA1 pyramidal neurons. Simplest observable tissue pattern: patches of dense reactive astrocytes sitting against local pyramidal-neuron fields in CA1.
In the donor table, median reactive hotspot exposure for the winning variation is 0.0011 and median astrocyte-hotspot exposure is 0.0035, consistent with the biomarker reading out the excess reactive component rather than absolute glial abundance.

## Next
Next local sweep: keep the same CA1 pyramidal-exposure setup but vary the exposure radius around hotspot cells (for example 50 px vs 90 px) while keeping the winning hotspot threshold fixed, to test whether the error pattern reflects an overly local or overly diffuse neuron-facing niche.
