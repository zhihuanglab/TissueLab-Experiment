## Summary
Tested the ca1_reactive_pyramidal_niche_atrophy family; candidate_variant_a won with selection_score=0.1178.

## Metrics
Winner candidate_variant_a (Reactive-adjacent radius = 60 px (about 30 µm), the baseline CA1 glial-neuron neighborhood scale.) had partial_r=-0.1178, selection_score=0.1178, loo_predictive_r=-0.0586, is_loo_gap=0.0591, penalty=0.0591, adjusted_score=0.0586.

Ranking of the other tested local variations:
- candidate_variant_b: partial_r=-0.0478, selection_score=0.0464, loo_predictive_r=-0.5185

## Findings
1. The best-performing signal came from CA1 pyramidal neurons that sit inside reactive-astrocyte neighborhoods: taking the median log-area shift relative to all CA1 pyramidal neurons isolates a niche-specific atrophy readout rather than global neuron loss.
2. The nearby alternative lost because tightening the radius changed adjacency counts and reduced the stability of the median without clearly enriching for a stronger injury state.
3. Error pattern:
- H20.33.018: outcome=-0.304, predicted=-0.097, feature=-0.0517 (overpredicted)
- H21.33.023: outcome=-0.155, predicted=-0.033, feature=-0.2784 (overpredicted)
- H21.33.035: outcome=-0.050, predicted=-0.161, feature=-0.0551 (underpredicted)
- H21.33.040: outcome=-0.033, predicted=-0.142, feature=-0.2942 (underpredicted)
- H20.33.037: outcome=-0.251, predicted=-0.143, feature=-0.1439 (overpredicted)

The most underpredicted donors were H21.33.035, H21.33.040, while the most overpredicted donors were H20.33.018, H21.33.023; both groups appear to deviate because CA1 reactive-adjacent soma shrinkage is only one component of their outcome residual.

## Rationale
This approach is biologically coherent because reactive astrocytes in CA1 already appear informative in the accepted panel, and soma shrinkage of neighboring pyramidal neurons is a plausible downstream injury phenotype. The winning radius likely balances perisomatic specificity against enough adjacent neurons for a stable donor-level median.

Likely additivity beyond the current panel: this candidate stays in the established CA1 reactive niche, but switches from composition/proximity to a morphology-severity readout, so it plausibly contributes partly new information if the evaluator confirms panel gain.

## Interpretation
The signal reflects CA1 Pyramidal Neuron morphology within a Reactive Astrocyte niche. The donor scalar is the reactive-adjacent minus global CA1 median log soma area. The simplest visible pattern would be smaller pyramidal neuron bodies preferentially where reactive astrocytes cluster in CA1.

Population: CA1 Pyramidal Neuron.  
Niche: CA1 Reactive Astrocyte-adjacent neighborhood.  
Feature summary: median(log area adjacent) minus median(log area of all CA1 pyramidal neurons).  
Observable pattern: preferential pyramidal soma shrinkage inside reactive astrocyte clusters.

## Next
Next, keep the same CA1 reactive-neighborhood injury theme but test whether a robust lower-tail summary (for example 25th percentile log area shift) outperforms the median, since the current errors suggest only a subset of adjacent neurons may be strongly damaged.
