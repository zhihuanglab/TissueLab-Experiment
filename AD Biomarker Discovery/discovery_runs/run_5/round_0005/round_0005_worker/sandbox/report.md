## Summary
Tested the CA1 peripyramidal astrocyte-lineage reactivity-index family, and the winning local variation was the 30 µm niche radius with a local winner selection score of **0.4866**.

## Metrics
Winning variation: `candidate_variant_a` (`ca1_peripyramidal_astro_reactivity_index_30um`).

- IS partial r: **-0.4866**
- Selection score: **0.4866**
- LOO predictive r: **0.3388**
- IS-LOO gap: **0.8255** with penalty **0.8255**
- Adjusted score: **-0.3388**

Additional diagnostics for the winner:
- n analyzable / total: **35 / 35**
- p-value: **0.0030**
- bootstrap median partial r: **-0.4925**
- bootstrap sign consistency: **0.9975**
- LOO unstable donors: **1** (max shift **0.1015**)

Ranking of tested nearby variations:
1. `candidate_variant_a` — CA1 peripyramidal astrocyte-lineage reactivity index using a 30 um radius around CA1 pyramidal neurons. (partial r **-0.4866**, selection **0.4866**, LOO **0.3388**)
2. `candidate_variant_b` — CA1 peripyramidal astrocyte-lineage reactivity index using a tighter 20 um radius around CA1 pyramidal neurons. (partial r **-0.4508**, selection **0.4508**, LOO **0.2861**)

## Findings
1. What worked and why
   - Restricting to **Astrocyte + Reactive Astrocyte** cells within the **CA1 pyramidal-neuron niche** worked better than the tighter nearby radius because it preserves the biology that has been strong across earlier rounds while changing the donor scalar from burden or size to a **state-balance fraction**.
   - The winning donor scalar is `ca1_peripyramidal_astro_reactivity_index_30um`: **reactive niche count / (reactive niche count + non-reactive astrocyte niche count)** among CA1 astrocyte-lineage cells whose centroids are within **30 µm** of any CA1 pyramidal neuron.
   - This likely captures a donor-level shift from homeostatic to reactive astrocyte state in the peripyramidal niche, which is biologically coherent with local neuronal stress.

2. What failed and why
   - The 20 µm variant was directionally similar but weaker, suggesting the signal is not confined to only the immediately juxtaneuronal ring; it benefits from a slightly broader peripyramidal neighborhood.
   - The family still targets the same CA1 reactive-astrocyte process as the accepted panel, so some redundancy is likely. The ratio summary helps, but it may still overlap with earlier burden and crowding readouts.
   - A too-tight niche probably drops informative lineage cells that still belong to the same local glial response field, reducing stability.

3. Error pattern: which donors are consistently wrong and what they share
   - The most unstable leave-one-out donors were: **H21.33.006**.
   - Largest underpredictions: **H21.33.006, H20.33.038, H20.33.018**.
   - Largest overpredictions: **H21.33.018, H21.33.035, H21.33.033**.
   - Overall pattern: some donors with severe decline still look less reactive by this ratio alone, while others have high local astrocyte reactivity balance without equally severe memory decline. That suggests the index tracks a real local glial state, but not the entire downstream severity spectrum.

## Rationale
The best variation is biologically coherent because it keeps the same disease-relevant population and niche already enriched by earlier rounds — **astrocyte-lineage cells around CA1 pyramidal neurons** — and asks a sharper question: within that niche, how much of the astrocyte pool has shifted into the reactive state? It beat the 20 µm alternative because the peripyramidal response appears to extend beyond only the closest cell-cell contacts and is better summarized over a slightly broader 30 µm field.

Relative to the current panel, this ratio is plausibly somewhat additive because it measures **state composition** rather than pure burden, exposure, crowding, or cell size. But it still sits on the same CA1 astrocyte axis, so later panel evaluation should check whether it adds new information or mostly re-expresses existing reactive-astrocyte load.

## Interpretation
Biologically, the signal seems to mean: donors with worse memory decline tend to show a **higher reactive share among astrocyte-lineage cells in the CA1 pyramidal-neuron neighborhood**.

- Population: **CA1 Astrocyte + Reactive Astrocyte lineage**
- Niche: **within 30 µm of a CA1 Pyramidal Neuron centroid**
- Feature summary: **reactive fraction among peripyramidal astrocyte-lineage cells**
- Simplest observable pattern: **around the CA1 pyramidal layer, a larger share of nearby astrocytes appears shifted into a reactive morphology/state rather than a more baseline astrocyte state**

## Next
Stay in the same niche and test a **donor-normalized excess-reactivity** sweep next: for example, compare the peripyramidal reactive fraction against each donor's broader CA1 astrocyte-lineage reactive fraction, to see whether the predictive signal comes from specifically pyramidal-adjacent state enrichment rather than general CA1 astrocyte reactivity.
