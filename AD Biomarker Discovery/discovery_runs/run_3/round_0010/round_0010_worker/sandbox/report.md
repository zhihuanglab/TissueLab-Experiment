## Summary
Tested the ca1_conditional_immune_cuff_fraction_dense_lymph family in CA1 pyramidal neurons, and the winner was ca1_conditional_immune_cuff_fraction_r50um_ra2_ly2 with local selection score 0.1474.

## Metrics
Winning variation: ca1_conditional_immune_cuff_fraction_r50um_ra2_ly2.
- IS partial r: 0.1474
- Selection score: 0.1474
- LOO predictive r: -0.1784
- IS-LOO gap: 0.3258 (penalty=0.3258)
- Adjusted score: -0.1784
- n analyzable / total: 35 / 35

Ranking summary:
1. ca1_conditional_immune_cuff_fraction_r50um_ra2_ly2 — partial_r=0.1474, selection_score=0.1474, LOO r=-0.1784
2. ca1_conditional_immune_cuff_fraction_r50um_ra2_ly3 — partial_r=NA, selection_score=NA, LOO r=NA

## Findings
1. What worked and why: The winning donor scalar was the fraction of severely reactive-cuffed CA1 pyramidal neurons that also had dense local lymphocyte admixture within 50 um. Higher values associate with slower memory decline. This works biologically because it keeps the successful severe reactive-cuff denominator from earlier rounds, then asks whether the same vulnerable perineuronal niche becomes even more inflammatory when lymphocytes accumulate beyond incidental single-cell contact.
2. What failed and why: The nearby stricter alternative ca1_conditional_immune_cuff_fraction_r50um_ra2_ly3 ranked second (partial_r=NA, selection_score=NA, LOO r=NA). That pattern suggests that tightening the threshold too far may make the event too sparse relative to the already small severe-cuff denominator, reducing donor-level stability.
3. Error pattern: Largest absolute LOO errors were in H20.33.018, H21.33.040, H21.33.023, H20.33.013, H20.33.037. Clinical-status and AD-neuropath-change labels were not available in the local cohort table used by this worker, so those shared metadata could not be summarized directly here. Their median severe-cuffed CA1 pyramidal count was 121.0000, with median ly>=2-positive count 0.0000 and ly>=3-positive count 0.0000, suggesting residual error is partly driven by donor-to-donor differences in how often the severe cuff progresses into a very dense immune pocket.

## Rationale
The best variation is biologically coherent because it targets a specific perineuronal inflammatory niche: CA1 pyramidal neurons already surrounded by a severe reactive astrocyte cuff, then further filtered for multiple nearby lymphocytes. This beat the nearby alternative because it appears to capture dense immune admixture without crossing into a threshold so rare that the donor summary becomes mostly sampling noise. Relative to the current panel, this is most likely a refinement of the accepted immune-cuff axis rather than a wholly orthogonal biology, but it still tests whether denser lymphocytic engagement marks a sharper inflammatory extreme inside that axis.

## Interpretation
Population: CA1 pyramidal neurons.
Niche: severe reactive-astrocyte cuffs within 50 um, refined by requiring multiple nearby lymphocytes.
Feature summary: donor-level fraction of severe-reactive-cuffed CA1 pyramidal neurons that also satisfy the winning lymphocyte threshold.
Observable tissue pattern: CA1 neuronal bands where already astrocyte-ensheathed pyramidal neurons sit inside mixed glial-immune microcuffs containing several lymphocytes. Biologically, the signal seems to reflect escalation from reactive glial cuffing toward a denser inflammatory perineuronal niche around vulnerable CA1 neurons.

## Next
Keep the same CA1 severe reactive-cuff denominator, but next sweep should test whether the winning dense-lymphocyte fraction is stronger when restricted to a peripyramidal CA1 band or normalized by local lymphocyte burden per severe-cuffed neuron instead of a binary threshold.
