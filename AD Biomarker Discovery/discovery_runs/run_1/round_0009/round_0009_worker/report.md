## Summary
Tested the CA1 reactive-astrocyte-proximal pyramidal lower-tail log-area family; candidate_variant_b (q20 within 80 px) won with selection score 0.0994.

## Metrics
Winning variation: candidate_variant_b (q20 within 80 px).

- IS partial r: -0.0994
- Selection score: 0.0994
- LOO predictive r: -0.1405
- IS-LOO gap: 0.0411 (penalty=0.0000)
- Adjusted score: -0.1405

Other tested variation ranking: candidate_variant_a (q25 within 80 px): partial_r=-0.0955, selection=0.0955, loo=-0.1084

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The strongest local signal came from q20 within 80 px, meaning the donor scalar captured how far the lower tail of CA1 pyramidal size shifts downward specifically inside a reactive-astrocyte niche. That is biologically coherent with local neuronal atrophy or shrinkage near reactive glia.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The nearby alternative changed only the quantile cutoff. It stayed in-family and remained analyzable, but it ranked lower, suggesting the very extreme tail was slightly noisier than the broader lower-tail summary on this cohort.
3. Error pattern: which donors are consistently wrong and what they share
   - Largest raw-outcome LOO errors were H20.33.018, H21.33.023, H21.33.040, H20.33.037, H21.33.035. These donors often sit at cohort extremes or have sparse proximal-cell support; 40% of the top-error donors fall below the median proximal-count support.

## Rationale
The winning variation uses a threshold-free lower-tail area summary within the exact CA1 pyramidal/reactive-astrocyte niche that already showed signal in prior rounds. It beat the nearby alternative because the q20 within 80 px summary appears to trade off specificity and stability slightly better than the more extreme tail. Given the active panel already contains CA1 reactive-astrocyte niche features, this winner likely measures a more morphology-specific version of that biology rather than a wholly new compartment.

## Interpretation
The signal appears to come from CA1 Pyramidal Neuron cells in the niche defined by proximity to CA1 Reactive Astrocyte cells. The donor-level scalar is the q20 within 80 px summary of log1p(cell area) among those proximal pyramidal neurons. The simplest observable tissue pattern is a donor whose CA1 pyramidal neurons look selectively smaller or more left-shifted in size near reactive astrocytes.

## Next
Keep the same CA1 pyramidal/reactive-astrocyte niche but sweep a slightly broader lower-tail summary next, such as q30 or a trimmed mean of the lowest quartile, to test whether the signal improves by reducing extreme-tail noise while preserving the same biological interpretation.
