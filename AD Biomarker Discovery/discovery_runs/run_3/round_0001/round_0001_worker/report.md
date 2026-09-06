## Summary
One sentence: tested astrocytic reactivity fractions in CA1 versus a broader pyramidal-field pool; candidate_variant_a won with local selection score 0.4914.

## Metrics
Winning variation: candidate_variant_a (CA1) with partial r -0.4914, selection score 0.4914, LOO predictive r 0.3598, IS-LOO gap 0.1316, penalty 0.1316, adjusted score 0.3598.
Other tested variations: candidate_variant_b had selection score 0.4629, partial r -0.4629, and LOO r 0.3422.

## Findings
1. What worked and why: The CA1-only astrocytic reactivity fraction worked best, consistent with CA1 being a memory-critical and AD-vulnerable niche where reactive gliosis may track local tissue stress more specifically than a broader hippocampal pool.
2. What failed and why: Pooling CA1 with CA2 and CA3 diluted the signal, implying that the strongest association is not a generic pyramidal-field astrocytic response but a more CA1-focused shift in astrocyte state.
3. Error pattern: Largest absolute LOO errors were H20.33.018, H21.33.023, H21.33.006, H21.33.018, H20.33.046. 0% of these donors were labeled Dementia. 0% had Intermediate/High AD neuropathologic change. Their median Braak/CERAD numeric values were 4.0/3.0. This suggests the astrocytic fraction alone does not fully capture how severe pathology translates into memory decline for these donors.

## Rationale
CA1 is selectively vulnerable in AD-related hippocampal degeneration, so a within-lineage reactive-astrocyte fraction in that region is biologically coherent as a compact readout of local gliosis. It beat the broader pooled alternative because adding CA2/CA3 likely mixed in less vulnerable tissue and reduced anatomical specificity.
If relevant to the future panel, this feature appears interpretable and auditable, but whether it adds new information beyond later panel members will depend on whether the eventual panel contains other glial or hippocampal-region features.

## Interpretation
The signal most plausibly reflects astrocytic gliosis. Population: Astrocyte plus Reactive Astrocyte lineage. Niche: CA1. Feature summary: Reactive Astrocyte / (Astrocyte + Reactive Astrocyte) at the donor level. Simplest observable pattern: donors with worse memory decline tend to show a larger share of astrocytes classified as reactive within the selected hippocampal niche.

## Next
Next sweep: keep the astrocyte/reactive-astrocyte lineage and CA1 gate fixed, but compare CA1 reactive-astrocyte fraction against CA1 reactive-astrocyte density or raw reactive count to test whether high-error donors are missing burden information rather than state-composition information.
