## Summary
Tested the ca1_reactive_astrocyte_enrichment family in CA1; candidate_variant_a won with selection score 0.4914.

## Metrics
Winning variation: candidate_variant_a (CA1 Reactive Astrocyte / (Reactive Astrocyte + Astrocyte)).
- IS partial r: -0.4914
- Selection score: 0.4914
- LOO predictive r: 0.3598
- IS-LOO Gap: 0.1316
- Gap penalty: 0.0000
- Adjusted Score: 0.4914
- Coverage: 35/35

Other tested variations ranked by selection score: candidate_variant_b (0.4217).

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The winning signal was the CA1 reactive astrocyte enrichment scalar `ca1_reactive_astrocyte_enrichment__reactive_over_astroglial`. It captures how much of the local CA1 astroglial compartment is shifted toward the reactive astrocyte state, which is biologically coherent with gliosis localized to the hippocampal subfield already implicated by the accepted CA1 pyramidal burden marker.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The losing variant diluted the signal by using all classified CA1 cells in the denominator. That mixes reactive astrocytes with large donor-to-donor shifts in neuronal and oligodendroglial abundance, so it is less specific to astroglial activation and more entangled with overall CA1 composition.
3. Error pattern: which donors are consistently wrong and what they share
   - Largest absolute LOO prediction errors occurred in donors H20.33.018, H21.33.023, H21.33.006. These donors likely share CA1 states where reactive gliosis does not move in lockstep with memory decline, suggesting heterogeneity between injury-response burden and downstream clinical effect.

## Rationale
The best variation is biologically coherent because it asks a state question inside the relevant lineage and region: within CA1 astroglia, what fraction is reactive rather than homeostatic astrocyte? That is a cleaner readout of local gliosis than counting reactive astrocytes against every CA1 cell class. It beat the nearby alternative because lineage normalization removes broad tissue-composition shifts and keeps the scalar anchored to astroglial state. Relative to the current panel's CA1 pyramidal burden feature, this candidate plausibly adds information about glial response rather than simple neuronal depletion.

## Interpretation
The signal appears to mean CA1-localized astroglial reactivity. Population: reactive astrocytes within the CA1 astroglial lineage. Niche: CA1 hippocampal annotation polygons. Feature summary: donor-level reactive astrocyte fraction within CA1 astroglia. Simplest observable tissue pattern: donors with higher biomarker values likely show CA1 fields where a larger share of astroglial cells are classified as reactive rather than resting astrocytes.

## Next
Run one nearby sweep that keeps CA1 and reactive astrocytes fixed but tests whether restricting the denominator to CA1 astroglia adjacent to pyramidal-neuron-rich neighborhoods improves the donors now dominating the LOO error pattern.
