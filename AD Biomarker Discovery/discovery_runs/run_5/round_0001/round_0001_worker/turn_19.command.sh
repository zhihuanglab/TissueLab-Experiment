cat > /scratch/report.md <<'MD'
## Summary
Reactive-astrocyte fraction among astrocyte-lineage cells was tested in a CA1-local versus CA1+CA2 gate; the CA1-only variation (`candidate_variant_a`) won with a local selection score of 0.4914.

## Metrics
Winning variation: `candidate_variant_a` (`ca1_reactive_astro_fraction`).

- IS partial r: -0.4914
- Selection score: 0.4914
- LOO predictive r: 0.3598
- IS-LOO Gap: 0.1316 (penalty=0.0000)
- Adjusted Score: 0.4914
- N analyzable: 35/35 donors

Ranking of tested variations:
1. `candidate_variant_a` — CA1 only, partial r = -0.4914, selection score = 0.4914, LOO predictive r = 0.3598
2. `candidate_variant_b` — CA1+CA2, partial r = -0.4801, selection score = 0.4801, LOO predictive r = 0.3507

The two variants are close, but adding CA2 slightly weakened both the in-sample partial signal and the LOO diagnostic, so the signal appears sharper in CA1 than in the broader pyramidal-field gate.

## Findings
1. **What worked and why**  
   The simple donor-level scalar worked: reactive astrocytes as a fraction of the astrocyte lineage inside CA1 tracked worse `slope_zmem0` after controlling for age, Braak, CERAD, and sex. The negative partial r means donors with a higher CA1 reactive-astrocyte share tended to have more negative memory decline slopes. This is biologically coherent for a memory-vulnerable hippocampal field, where local gliosis can mark tissue stress or injury not fully summarized by standard neuropathology staging.
2. **What failed and why**  
   Broadening the gate from CA1 to CA1+CA2 lost a small amount of signal. That suggests the reactive-astrocyte enrichment is not just a generic hippocampal astrocytosis measure; it seems more localized to CA1. The broader gate likely diluted CA1-specific variation by mixing in CA2 tissue with a somewhat different glial composition or disease relationship.
3. **Error pattern**  
   The largest underpredictions were donors `H20.33.018` and `H20.33.037`, both with high CA1 reactive-astrocyte fractions and substantially worse memory decline than predicted, suggesting additional damaging processes beyond this glial marker. The largest overprediction was `H21.33.006`, which had a very high CA1 reactive-astrocyte fraction (0.8106) but only mild decline, implying reactive astrocytosis can sometimes be disproportionate to memory loss. Donors `H21.33.023` and `H21.33.018` were also underpredicted despite relatively low fractions, pointing to non-astrocytic or non-CA1 contributors in those cases.

## Rationale
The winning approach is biologically coherent because it measures a lineage-normalized state shift rather than raw abundance. Using `Reactive Astrocyte / (Astrocyte + Reactive Astrocyte)` focuses on how much of the local astrocyte compartment appears reactive, which is more interpretable than total astrocyte count and less confounded by overall cellularity. CA1 beat CA1+CA2 because the scientific question was specifically about a memory-vulnerable field, and the data support that narrower niche: the broader adjacent region did not add signal and instead slightly diluted it.

Because the current accepted panel is empty, this looks like a plausible seed biomarker. It is simple, fully replayable from raw cell labels plus region annotations, and likely captures biology distinct from pure plaque/tangle stage because the association persists after Braak and CERAD adjustment.

## Interpretation
The signal seems to mean **CA1-local reactive gliosis intensity within the astrocyte lineage**.

- **Population:** Astrocyte-lineage cells (`Astrocyte` and `Reactive Astrocyte`)
- **Niche:** CA1 hippocampal annotation gate
- **Feature summary:** donor-level fraction `Reactive Astrocyte / (Astrocyte + Reactive Astrocyte)`
- **Simplest observable tissue pattern:** donors with more CA1 astrocytes shifted into the reactive state tend to show worse memory decline

Operationally, this corresponds to slides where the CA1 field contains a larger share of cells classified as reactive astrocytes relative to quiescent astrocytes, consistent with local injury-response gliosis.

## Next
Run the next local sweep within the same glial-state family by keeping CA1 fixed and changing only the donor summary: compare the current lineage fraction against a CA1 reactive-astrocyte density per CA1 area/cellularity-normalized denominator, because the current errors suggest some donors have unusually high reactive-state fraction without proportionally severe decline.
MD
test -s /scratch/report.md
