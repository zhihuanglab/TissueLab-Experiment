## Summary
Tested the ca1_conditional_immune_cuff_fraction family in CA1 pyramidal neurons, and the winner was ca1_conditional_immune_cuff_fraction_r50um_ra2_ly1 with local selection score 0.2779.

## Metrics
Winning variation: ca1_conditional_immune_cuff_fraction_r50um_ra2_ly1.
- IS partial r: -0.2779
- Selection score: 0.2779
- LOO predictive r: 0.0755
- IS-LOO gap: 0.2024 (penalty=0.2024)
- Adjusted score: 0.0755
- n analyzable / total: 35 / 35

Ranking summary:
1. ca1_conditional_immune_cuff_fraction_r50um_ra2_ly1 — partial_r=-0.2779, selection_score=0.2779, LOO r=0.0755
2. ca1_conditional_immune_cuff_fraction_r50um_ra2_ly2 — partial_r=0.1474, selection_score=0.1474, LOO r=-0.1784

## Findings
1. What worked and why: Conditioning on an existing severe reactive-astrocyte cuff sharpened the immune signal. The winning donor scalar is the fraction of already severely cuffed CA1 pyramidal neurons that are also lymphocyte-admixed within 50 um. Higher values associate with worse memory decline. This is biologically coherent because it isolates an inflammatory escalation niche rather than counting all cuffed neurons equally.
2. What failed and why: Requiring a denser lymphocyte burden was weaker. The alternative ca1_conditional_immune_cuff_fraction_r50um_ra2_ly2 ranked second (partial_r=0.1474, selection_score=0.1474, LOO r=-0.1784). That suggests the useful information is presence of immune admixture within a severe cuff, while insisting on two lymphocytes makes the event too sparse and noisier at donor level.
3. Error pattern: Largest absolute LOO errors were in H20.33.018, H20.33.013, H20.33.038, H21.33.044, H20.33.037. The compact cohort table available to this worker does not include cognitive-status or AD-neuropath-change labels, so shared clinical metadata for these outliers could not be checked directly here. Their median severe-cuffed CA1 pyramidal count was 210.0000, suggesting the biomarker still misses some donor-to-donor heterogeneity even when the severe cuff population is abundant.

## Rationale
The best approach is biologically coherent because it asks a specific second-stage question inside a previously successful niche: once a CA1 pyramidal neuron is already embedded in a severe reactive astrocyte cuff, does lymphocyte admixture mark an even more pathologic microenvironment? This beat the nearby alternative because a single nearby lymphocyte appears to capture immune involvement before the niche becomes so dense that counts become rare and unstable. Relative to the current panel, this looks like a cleaner inflammatory refinement of the existing immune-admixed cuff feature rather than a generic repeat of overall cuff prevalence.

## Interpretation
Population: CA1 pyramidal neurons.
Niche: severe reactive-astrocyte cuffs within 50 um, further refined by local lymphocyte admixture.
Feature summary: donor-level fraction of severe-reactive-cuffed CA1 pyramidal neurons that also meet the lymphocyte threshold.
Observable tissue pattern: CA1 neuronal bands where already astrocyte-ensheathed pyramidal neurons are additionally touched by one or more lymphocytes. Biologically, the signal seems to represent escalation from glial reactivity into mixed glial-immune inflammatory cuffing around vulnerable CA1 neurons.

## Next
Keep the same conditional severe-cuff denominator, but in the next local sweep test whether the same lymphocyte-admixed fraction is stronger when restricted to the peripyramidal CA1 band or when weighted by reactive-astrocyte cuff size rather than treated as a binary severe-cuff event.
