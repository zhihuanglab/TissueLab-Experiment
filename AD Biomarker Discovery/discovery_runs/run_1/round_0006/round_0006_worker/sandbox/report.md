## Summary
Tested the CA1 reactive-astrocyte-proximal pyramidal soma-area family; the 80 px radius median-log-area variant won locally with selection score 0.0423.

## Metrics
Winning variation: `median_log_area_80px` (`ca1_pyramidal_reactive_astro_proximal_median_log_area_80px`).

- IS partial r: -0.0423
- Selection score: 0.0423
- LOO predictive r: -0.2999
- IS-LOO Gap: 0.0000
- Penalty: 0.0000
- Adjusted Score: 0.0423

Other tested variation:
- `median_log_area_60px`: partial r -0.0033, selection score 0.0033, LOO predictive r -0.5948

Coverage was complete for both variants (35/35 donors analyzable). For the winning 80 px niche, the selected-cell count was not sparse: median 419 CA1 pyramidal neurons per donor were within the reactive-astrocyte gate, with range 54-1118.

## Findings
1. What worked and why  
   The only thing that worked at all was broadening the reactive-astrocyte proximity gate from 60 px to 80 px. That slightly stabilized the donor summary by averaging morphology over a larger CA1 pyramidal niche, suggesting the niche definition itself is feasible and the donor scalar is measurable with good coverage.
2. What failed and why  
   The core hypothesis failed: donor-level median log soma area of CA1 pyramidal neurons near reactive astrocytes carried essentially no confound-adjusted in-sample association with `slope_zmem0` (partial r only -0.0423). This implies that in this dataset, reactive-astrocyte-adjacent CA1 pyramidal neurons do not show a strong monotonic area shift that tracks memory decline, even though the same niche was informative in prior composition/coverage rounds.
3. Error pattern: which donors are consistently wrong and what they share  
   Largest LOO errors came from H20.33.018, H21.33.023, H21.33.040, H20.33.037, and H20.33.038. These donors are mostly moderate-to-high Braak cases (often 4-6), so the model tends to miss donors with heavier pathology burden, but the errors do not collapse into one sex or one obvious feature-value extreme. This looks more like morphology/outcome decoupling than a simple subgroup coverage issue.

## Rationale
The best variation is still biologically coherent: it keeps the previously supported CA1 pyramidal/reactive-astrocyte niche and asks whether a directly degenerative morphology readout, pyramidal contour area, adds information. The 80 px version beat 60 px because it likely reduced donor-level noise by including more niche-proximal cells while staying local to the same astrocyte-conditioned microenvironment. Even so, the effect size was tiny, so this morphology summary looks unlikely to add meaningfully beyond the current panel.

## Interpretation
Biologically, this signal appears to mean very little for `slope_zmem0` in the training cohort.  
Population: CA1 pyramidal neurons.  
Niche: cells lying within 80 px of CA1 reactive astrocytes.  
Feature summary: donor-level median `log1p(contour area_px^2)`.  
Simplest observable tissue pattern: a donor whose CA1 pyramidal somata near reactive astrocytes are slightly smaller would score higher in the negative direction, but the association is so weak that no reliable visual pattern is supported here.

## Next
Stay on the same CA1 pyramidal/reactive-astrocyte niche but switch to a nearby morphology variation that is more sensitive to shape distortion than bulk size, such as elongation, circularity, or donor tail-fraction of abnormally small pyramidal contours within the 80 px gate.
