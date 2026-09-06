## Summary
One sentence: tested the `ca1_reactive_hotspot_lymphocyte_enrichment` family, the winning variation was `reactive_minus_astro_hotspot_local_lymphocyte_fraction_r70px_k3`, and the local winner selection score was `0.0916`.

## Metrics
Winning variation `reactive_minus_astro_hotspot_local_lymphocyte_fraction_r70px_k3` (`ca1_reactive_hotspot_lymphocyte_enrichment__reactive_minus_astro_hotspot_local_lymphocyte_fraction_r70px_k3`): IS partial r = 0.1145, selection score = 0.0916, LOO predictive r = 0.0486, IS-LOO gap = 0.0659 (penalty=0.0000), adjusted score = 0.0486.
Other ranked variations:
- `reactive_minus_astro_hotspot_local_lymphocyte_fraction_r90px_k3`: partial_r=0.0874, selection_score=0.0699, loo_predictive_r=0.0532

## Findings
1. What worked and why (tie to the biological meaning of the target): The tighter 70 px readout preserved the most focal immune niche contrast around reactive hotspots. Biologically, the feature asks whether CA1 reactive astrocyte microclusters sit in a more lymphocyte-rich microenvironment than matched homeostatic astrocyte microclusters.
2. What failed and why (specific to the chosen hypothesis and what went wrong): The broader reactive_minus_astro_hotspot_local_lymphocyte_fraction_r90px_k3 window likely diluted a sparse lymphocyte signal by averaging in surrounding CA1 background. The losing local sweep stayed in the same niche family but changed only the measurement radius.
3. Error pattern: The largest LOO errors were H20.33.030, H20.33.018, H20.33.037; several were more impaired than predicted, suggesting decline not fully captured by focal lymphocyte enrichment; others carried the niche signal without equally severe decline, suggesting inflammatory-context signal can decouple from clinical slope.

## Rationale
The best variation is biologically coherent because it isolates a specific inflammatory niche: hotspot-positive Reactive Astrocytes in CA1, benchmarked against hotspot-positive Astrocytes from the same region, and summarized as a donor-level delta in local lymphocyte fraction. That comparison controls away global cell-density and region-size effects better than a whole-CA1 lymphocyte summary. It beat the nearby alternative because its radius gave the better balance between focal specificity and stability for a rare immune-cell context. It may add information beyond the current panel because prior kept features in this hotspot series emphasized neuronal loss, oligodendrocyte depletion, and corpora-amylacea context rather than explicit lymphocyte enrichment.

## Interpretation
The signal seems to mean that some donors have CA1 reactive astrocyte hotspots embedded in a more immune-infiltrated local niche than ordinary astrocyte hotspots. Population: Reactive Astrocyte versus Astrocyte hotspot-positive cells in CA1. Niche: local CA1 neighborhoods within the winning radius around those hotspot centers. Feature summary: donor-level mean local lymphocyte fraction around reactive hotspots minus the same mean around astrocyte hotspots. Simplest observable pattern: focal reactive astrocyte microclusters that appear to sit nearer small lymphocyte pockets than nearby baseline astrocyte clusters.

## Next
Keep the same CA1 reactive-hotspot framework, but next sweep whether the immune context is better captured by a tighter cell-state gate inside the same family, for example hotspot-local lymphocyte presence versus fraction or a reactive-hotspot subset with stronger same-class crowding.
