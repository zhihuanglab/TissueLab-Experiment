## Summary
Tested the CA1 peripyramidal reactive-astrocyte/lymphocyte contact family; the tighter 20 um variation (`candidate_variant_b`) won locally with selection score 0.1762.

## Metrics
Winning variation: `candidate_variant_b` (`ca1_peripyramidal_reactive_astro_lymphocyte_contact_fraction_20um`).

- IS partial r: 0.1762
- Selection score: 0.1762
- LOO predictive r: 0.0171
- IS-LOO Gap: 0.1591 (penalty=0.0046)
- Adjusted Score: 0.1717
- Coverage: 35/35 donors

Ranking of tested variations:
1. `candidate_variant_b` — partial_r 0.1762, selection_score 0.1762, loo_predictive_r 0.0171, adjusted_score 0.1717
2. `candidate_variant_a` — partial_r 0.1012, selection_score 0.1012, loo_predictive_r -0.0461, adjusted_score 0.1012

So the more local 20 um lymphocyte radius beat the baseline 30 um radius, while the broader radius diluted the already weak single-feature signal.

## Findings
1. **What worked and why**
   - Restricting the signal to **CA1 Reactive Astrocytes** that are already **peripyramidal** sharpened the biology to the validated astrocyte-neuron niche rather than global inflammation.
   - Within that niche, the **20 um lymphocyte contact fraction** performed better than 30 um, suggesting the informative pattern is a very local inflammatory adjacency rather than a broad lymphocyte-rich halo.
   - The signal appears to come from a small subset of donors with nonzero contact scores, especially donors such as `H20.33.013`, `H20.33.030`, `H21.33.043`, and `H21.33.044`, where some neuron-adjacent reactive astrocytes also sit close to CA1 lymphocytes.

2. **What failed and why**
   - The whole family is extremely sparse. The winning 20 um feature was nonzero in only 5 of 35 donors, so most donors collapse to the same value of 0.
   - Because of that sparsity, the feature has almost no out-of-sample usefulness despite a small positive in-sample partial correlation: LOO predictive r was only 0.0171.
   - The 30 um variation likely admitted more diffuse, less specific lymphocyte proximity and therefore weakened the local inflammatory-niche definition instead of improving coverage.

3. **Error pattern**
   - The largest misses were `H20.33.018`, `H21.33.040`, `H21.33.023`, `H20.33.037`, `H21.33.035`, and `H20.33.038`.
   - These donors mostly share the same pattern: **feature = 0 despite many peripyramidal reactive astrocytes**. For example, `H20.33.018` had 424 peripyramidal reactive astrocytes and 10 CA1 lymphocytes but still a 20 um contact fraction of 0.0.
   - The model therefore underestimates donors with worse memory decline when they lack this very rare lymphocyte-contact event. This suggests the inflammatory niche is real in a few cases but too infrequent to track the broader decline signal.

## Rationale
The best variation is biologically coherent because it asks a sharper question than the existing astrocyte-only panel: not just whether CA1 reactive astrocytes cluster around pyramidal neurons, but whether those neuron-adjacent reactive astrocytes also sit in a **local lymphocyte-associated inflammatory microenvironment**. The 20 um version beat 30 um because true immune-cell adjacency is likely highly focal; expanding the radius to 30 um appears to add nonspecific nearby lymphocytes without improving generalization.

That said, the diagnostic profile suggests this feature is unlikely to add much new information beyond the current CA1 astrocyte panel. Its in-sample gain is small and its LOO signal is nearly null, which is what you expect from a biologically plausible but very low-prevalence event.

## Interpretation
Biologically, the signal seems to represent a **rare inflammatory overlay on the established CA1 reactive astrocyte–pyramidal niche**.

- **Population:** CA1 Reactive Astrocytes
- **Niche:** peripyramidal subset defined by a CA1 Pyramidal Neuron within 30 um
- **Feature summary:** fraction of those niche-defined reactive astrocytes with at least one CA1 Lymphocyte within 20 um
- **Simplest observable pattern:** a small subset of reactive astrocytes hugging CA1 pyramidal neurons also have a nearby lymphocyte in the same local microenvironment

When present, this may correspond to focal immune engagement around already stressed CA1 neuron-astrocyte interfaces, but in this cohort it is too rare to serve as a robust standalone donor biomarker.

## Next
Stay in the same CA1 peripyramidal reactive-astrocyte/lymphocyte branch, but replace the sparse binary contact fraction with a **continuous nearest-lymphocyte distance summary** (for example the median or lower-tail distance from peripyramidal reactive astrocytes to the nearest CA1 lymphocyte), because 20 um beat 30 um yet both thresholded versions fail mainly by assigning 0 to most donors.
