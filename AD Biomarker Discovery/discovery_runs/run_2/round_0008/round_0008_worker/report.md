## Summary
One sentence: tested the ca1_reactive_hotspot_oligodendrocyte_depletion family, candidate_variant_b won, and its local winner selection score was 0.3976.

## Metrics
Winning variation: candidate_variant_b (`ca1_reactive_hotspot_oligodendrocyte_depletion__candidate_variant_b`).
- IS partial r: 0.3976
- Selection score: 0.3976
- LOO predictive r: 0.2580
- IS-LOO Gap: 0.1396
- Penalty: 0.0000
- Adjusted Score: 0.2580
- Coverage: 1.0000 (35/35)

Other ranked variations: candidate_variant_a partial_r=0.3296, selection=0.3296, loo_r=0.2348

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The family remained biologically coherent because it stayed inside CA1 reactive-versus-homeostatic astroglial hotspots and converted that niche contrast into a simple donor scalar: median CA1 oligodendrocyte neighbor count around reactive centers minus the same summary around astrocyte centers.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The tested depletion contrast produced moderate rather than dominant local single-feature signal, so replacing pyramidal-neuron readout with oligodendrocyte neighbor counts recovered a real niche effect but not an obviously overwhelming one. That suggests the hotspot frame may capture a myelin-support dimension that is present yet still partly intertwined with the broader CA1 injury program.
3. Error pattern: which donors are consistently wrong and what they share
   - Largest LOO misses: H20.33.018, H21.33.023, H20.33.038, H21.33.001, H21.33.040. Most large-error donors used hotspot fallback on at least one astroglial side.

## Rationale
The winning variation beat the nearby alternative because its hotspot threshold produced the better local single-feature selection score while preserving donor coverage. Biologically, the approach asks whether CA1 reactive astrocyte micro-hotspots sit in oligodendrocyte-poor niches relative to matched homeostatic astrocyte hotspots. Using the median local oligodendrocyte count is auditable and robust to outliers, and the stronger k=4 hotspot threshold slightly sharpened the donor-level effect versus k=3 while keeping full coverage. The local winner is strong enough that it may contribute information beyond the current CA1 gliosis/neuron-loss panel.

## Interpretation
This signal is defined in the CA1 astroglial hotspot niche. The population is Reactive Astrocyte versus Astrocyte centers, the niche is a 70 px CA1 hotspot neighborhood, the donor-level feature summary is reactive-center median oligodendrocyte neighbors minus astrocyte-center median oligodendrocyte neighbors, and the simplest observable tissue pattern is whether reactive astrocyte clusters are surrounded by visibly fewer nearby oligodendrocytes than matched homeostatic astrocyte clusters. In this run, more positive values of the biomarker correspond to stronger relative oligodendrocyte depletion around reactive hotspots.

## Next
One specific suggestion for the next local sweep based on the error pattern and which nearby variations won or lost: keep the same CA1 reactive-hotspot framework but test a ratio-normalized oligodendrocyte exposure summary (for example reactive/astro median local oligodendrocyte count or oligo fraction among all neighbors) to determine whether the absolute-count version lost signal mainly to donor-level cellularity differences.
