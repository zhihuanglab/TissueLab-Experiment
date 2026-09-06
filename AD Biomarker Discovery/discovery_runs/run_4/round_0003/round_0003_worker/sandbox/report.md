## Summary
Tested the `ca1_pyramidal_reactive_astrocyte_neighbor_fraction` family in CA1; `r6_neighbor_fraction` was the local winner with selection score 0.1986, but the signal was sparse and its LOO predictive r was only 0.0102.

## Metrics
Winning variation: `r6_neighbor_fraction`

- IS partial r: 0.1986
- Selection score: 0.1986
- LOO predictive r: 0.0102
- IS-LOO gap: 0.1884
- Penalty: 0.0192
- Adjusted score: 0.1794

Other tested variation:

- `r4_neighbor_fraction`: partial_r=0.0823, selection_score=0.0823, loo_predictive_r=-0.2714, adjusted_score=0.0627

The ranking favored radius 6 over radius 4 because the broader niche produced the larger confound-adjusted in-sample signal, even though neither radius generalized well in LOO.

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The family is biologically targeted: it asks whether the existing CA1 pyramidal-loss and CA1 reactive-astrocyte signals become sharper when they are forced to co-occur in the same micro-niche.
   - The radius-6 variant did slightly better because it captured a few more nonzero donor events at this coordinate scale: 4 donors were nonzero at radius 6 versus 3 at radius 4.
   - The donor-level scalar is still interpretable: fraction of CA1 pyramidal neurons with at least one nearby CA1 reactive astrocyte.

2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The dominant failure was sparsity. For the winning radius-6 feature, 31/35 donors were exactly 0 and the maximum donor value was only 0.000647.
   - Radius 4 was even sparser and produced a negative LOO predictive correlation, so tightening the niche further made the feature too close to a near-constant binary rarity signal.
   - This suggests the planned radii are too small relative to observed CA1 centroid spacing in this dataset, so the niche definition is biologically plausible but operationally too strict.

3. Error pattern: which donors are consistently wrong and what they share
   - The largest raw LOO misses were H20.33.018, H21.33.040, H21.33.023, H20.33.037, and H21.33.001.
   - Four of those five donors had feature value exactly 0 despite substantial CA1 pyramidal and reactive-astrocyte counts, so the model defaulted toward confound-driven mean predictions and missed donors with unusually steep or unusually mild decline.
   - The common pattern is not missing CA1 tissue or missing target populations; it is that the current niche rule almost never fires even when both populations are abundant in CA1.

## Rationale
The best variation remains biologically coherent because a neuron-glia injury niche is a more mechanistic readout than either CA1 pyramidal fraction or CA1 reactive astrocyte enrichment alone. Radius 6 beat radius 4 because it relaxed the neighborhood just enough to produce a slightly less degenerate donor distribution while staying within the same hypothesized local niche. Even so, the near-zero LOO predictive r suggests this feature is unlikely to add useful new information beyond the current accepted panel unless the neighborhood definition is broadened or re-parameterized.

## Interpretation
This signal is meant to represent a CA1-localized reactive gliosis niche around vulnerable pyramidal neurons.

- Population: CA1 pyramidal neurons
- Niche: immediate CA1 reactive-astrocyte neighborhood
- Feature summary: fraction of CA1 pyramidal neurons with at least one CA1 reactive astrocyte within the winning radius
- Simplest observable pattern: reactive astrocytes sitting directly adjacent to pyramidal neurons in CA1

In this run, the biology is plausible, but the measured pattern is so rare at radii 4-6 that the resulting donor scalar behaves more like a sparse event detector than a stable gradient of tissue injury.

## Next
Keep the same CA1 pyramidal–reactive astrocyte niche concept, but in the next local sweep test a materially larger neighborhood scale or a continuous nearest-distance summary, because the current binary radii of 4 and 6 coordinate units are too sparse to exploit the underlying biology.
