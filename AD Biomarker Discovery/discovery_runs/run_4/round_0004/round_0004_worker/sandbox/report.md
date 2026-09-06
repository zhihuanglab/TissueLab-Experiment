## Summary
Tested the ca1_reactive_astrocyte_hypertrophy family in CA1; median_log_area_shift won with selection score 0.3588.

## Metrics
Winner `median_log_area_shift` (`ca1_reactive_astrocyte_hypertrophy_median_log_area_shift`):
- IS partial r: -0.3588
- Selection score: 0.3588
- LOO predictive r: 0.2220
- IS-LOO Gap: 0.1368 (penalty=0.0000)
- Adjusted Score: 0.2220
- n_analyzable: 35/35

Other tested variation ranking:
- upper_tail_log_area_shift selection=0.3131, partial_r=-0.3131, loo_r=0.1599

## Findings
1. What worked and why: The median-shift version won, which suggests the relevant signal is a broad donor-level enlargement of CA1 reactive astrocytes relative to baseline astrocytes. The donor-level scalar is the within-CA1 log-area shift between Reactive Astrocytes and Astrocytes, so it tracks activation severity rather than raw abundance.
2. What failed and why: The upper-tail version appears too tail-sensitive on this cohort size, so a few extreme contours likely added noise instead of sharpening the donor summary.
3. Error pattern: largest held-out errors were H20.33.018 (obs=-0.304, pred=-0.123, feature=0.143, astro=406, reactive=1161); H21.33.023 (obs=-0.155, pred=-0.028, feature=-0.199, astro=932, reactive=449); H20.33.038 (obs=-0.237, pred=-0.130, feature=-0.251, astro=241, reactive=747). The biggest misses are not primarily missing-data cases; they look more like donors where reactive astrocyte size and memory decline are partly decoupled.

## Rationale
A median shift is biologically coherent when reactive enlargement is a widespread state change within the CA1 astroglial compartment, making it more stable than a tail-only readout. Restricting the comparison to CA1 sharpens the feature to the same hippocampal field where prior rounds already found pyramidal and reactive-astrocyte signal. Using CA1 Astrocytes as the within-slide reference means the biomarker asks whether reactive astrocytes are unusually enlarged beyond the baseline astrocyte size context of that same donor. That makes it more morphology-specific than abundance or proximity alone. Relative to the current panel, this feature looks biologically adjacent to the accepted CA1 reactive-astrocyte enrichment axis, but it is mechanistically more state-focused, so it still seems plausible as a source of partially nonredundant information.

## Interpretation
Biologically, the signal appears to mean stronger CA1 astrocyte activation marked by reactive-cell hypertrophy. The population is CA1 Reactive Astrocytes referenced against CA1 Astrocytes; the niche is the CA1 field itself; the donor summary is `median_log_area_shift` on log contour area; and the simplest observable tissue pattern is a donor having visibly larger reactive astrocyte contours in CA1 than its non-reactive astrocyte baseline.

## Next
Stay on the same CA1 reactive-astrocyte morphology axis and test one nearby sweep that changes only the donor aggregation: keep the winning CA1 reactive-vs-astrocyte size contrast, but compare this median_log_area_shift statistic with a robust trimmed-upper-tail summary to see whether the winner is capturing a broad state shift or a sparse hypertrophic subpopulation.
