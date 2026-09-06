## Summary
One sentence: tested the ca1_reactive_astrocyte_proximal_small_pyramidal_fraction family in CA1, candidate_variant_b won, and its local winner selection score was 0.3287.

## Metrics
Winner: candidate_variant_b with IS partial r 0.3287, selection score 0.3287, LOO predictive r 0.2259, IS-LOO gap 0.1028 (penalty 0.0000), adjusted score 0.2259. Other tested variations: candidate_variant_a (selection 0.2568).

## Findings
1. What worked and why (tie to the biological meaning of the target): The winning feature worked by summarizing the fraction of CA1 pyramidal neurons that are both reactive-astrocyte-proximal and in the donor-specific small-area tail. That is biologically coherent with a local degenerative/atrophic niche around reactive astrocytes rather than a global CA1 size shift.
2. What failed and why (specific to the chosen hypothesis and what went wrong): candidate_variant_a was weaker than the winner (selection 0.2568 vs 0.3287), consistent with its stricter tail cutoff being too sparse and missing moderately shrunken proximal pyramidal neurons.
3. Error pattern: The largest raw LOO errors were H20.33.018, H21.33.023, H20.33.037; these donors span cognitive states [unknown] with Braak [4, 5, 6] and CERAD [0, 3], suggesting the niche signal misses some pathology severity heterogeneity.

## Rationale
The best variation keeps the previously productive CA1 pyramidal–reactive astrocyte niche but converts the continuous area shift into a donor-normalized tail burden. Using the 33rd percentile beat the stricter 25th percentile cutoff, which suggests the signal is not confined to the most extreme atrophic tail and is better captured by a somewhat broader burden of shrunken proximal cells. Because the accepted panel already contains CA1 pyramidal abundance, reactive astrocyte lineage burden, and a proximal median-area term, this feature is most plausible as a thresholded tail version of the same niche morphology rather than a wholly orthogonal biology.

## Interpretation
The signal seems to mean that donors with higher burden of unusually small CA1 pyramidal neurons immediately adjacent to reactive astrocytes tend to have higher slope_zmem0. Population: CA1 pyramidal neurons. Niche: within 80 px of CA1 reactive astrocytes. Feature summary: fraction of proximal pyramidal neurons below a donor-specific CA1 pyramidal area quantile. Simplest observable pattern: patches of CA1 pyramidal neurons next to reactive astrocytes look shrunken relative to the donor's broader CA1 pyramidal population.

## Next
Run one more local sweep in the same CA1 pyramidal–reactive astrocyte niche, keeping the winning tail definition fixed while varying the proximity radius around 60–100 px to test whether the signal is truly immediate-neighbor atrophy or a broader peri-reactive-astrocyte field effect.
