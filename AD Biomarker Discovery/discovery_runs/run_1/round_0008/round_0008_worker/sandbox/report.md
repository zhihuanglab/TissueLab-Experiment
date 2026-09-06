## Summary
Tested a CA1 reactive-astrocyte-proximal small-pyramidal enrichment family; candidate_variant_b won with local selection score 0.2140.

## Metrics
Best variation: candidate_variant_b
- IS partial r: 0.2140
- Selection score: 0.2140
- LOO predictive r: 0.1864
- IS-LOO Gap: 0.0276 (penalty=0.0000)
- Adjusted Score: 0.2140
- Feature column: `ca1_small_pyramidal_enrichment_q20_80px`

Other tested variations:
- candidate_variant_a: partial_r=0.1963, selection_score=0.1963, loo_predictive_r=0.1673

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The winner worked best when the donor-level scalar asked whether morphologically small CA1 pyramidal neurons were selectively enriched inside the reactive-astrocyte-proximal compartment rather than globally. That contrast is aligned with a localized vulnerability signal rather than simple overall neuronal size or abundance.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The broader 25th-percentile threshold likely mixed the informative tail with less selective background size variation. It remained the same CA1 pyramidal/reactive-astrocyte niche, but the small-cell cutoff changed the bias-variance tradeoff inside a limited proximal compartment.
3. Error pattern: which donors are consistently wrong and what they share
   - Largest absolute LOO errors: H20.33.018 (outcome=-0.304, predicted=-0.103, ca1_small_pyramidal_enrichment_q20_80px=0.022); H21.33.023 (outcome=-0.155, predicted=-0.013, ca1_small_pyramidal_enrichment_q20_80px=0.058); H20.33.037 (outcome=-0.251, predicted=-0.130, ca1_small_pyramidal_enrichment_q20_80px=0.054); H21.33.040 (outcome=-0.033, predicted=-0.132, ca1_small_pyramidal_enrichment_q20_80px=0.060); H21.33.035 (outcome=-0.050, predicted=-0.148, ca1_small_pyramidal_enrichment_q20_80px=0.030). These errors tend to reflect donors where the niche contrast does not map cleanly onto the outcome despite measurable CA1 pyramidal and reactive astrocyte presence.

## Rationale
The winning variation is biologically coherent because it focuses on a specific population (CA1 pyramidal neurons), a specific niche (within 80 px of CA1 reactive astrocytes), and a donor-relative size-tail summary (20th-percentile within-donor small-cell threshold). It beat the nearby alternative because the chosen cutoff better isolated the donor-specific small-cell pattern without relying on an absolute area scale. Given the accepted panel already contains CA1 pyramidal fraction, reactive astrocyte lineage fraction, niche fraction, niche coverage, and proximal median log area, this enrichment contrast still looks plausibly additive because it measures selective concentration of the small-size tail between proximal and distal CA1 compartments rather than another absolute burden summary.

## Interpretation
The signal appears to represent localized enrichment of morphologically small, potentially stressed or atrophic CA1 pyramidal neurons around reactive astrocytes. Population: CA1 Pyramidal Neuron. Niche: CA1 pyramidal cells within 80 px of CA1 Reactive Astrocytes versus the distal CA1 background. Feature summary: proximal small-cell fraction minus distal small-cell fraction. Simplest observable pattern: donors with worse memory decline appear to show more shrinkage-like small pyramidal profiles concentrated near reactive astrocyte neighborhoods inside CA1.

## Next
Run the next local sweep around the same CA1 pyramidal/reactive-astrocyte niche but keep the winning proximal-minus-distal contrast and vary only the proximity radius (for example 60 px versus 100 px) to determine whether the signal is truly juxtacrine-local or a broader reactive-glial field effect.
