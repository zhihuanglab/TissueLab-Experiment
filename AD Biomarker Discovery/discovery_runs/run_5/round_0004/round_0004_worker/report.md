## Summary
Tested the CA1 peripyramidal reactive-astrocyte hypertrophy family, and the winning local variation was the 30 µm peripyramidal **median** reactive-astrocyte area with a local winner selection score of **0.3755**.

## Metrics
Winning variation: `candidate_variant_a` (`ca1_peripyramidal_reactive_astro_area_median_30um_um2`).

- IS partial r: **-0.3755**
- Selection score: **0.3755**
- LOO predictive r: **0.0995**
- IS-LOO gap: **0.4750** with penalty **0.4750**
- Adjusted score: **-0.0995**

Additional diagnostics for the winner:
- n analyzable / total: **35 / 35**
- p-value: **0.0262**
- bootstrap median partial r: **-0.3852**
- bootstrap sign consistency: **0.9225**
- LOO unstable donors: **1** (max shift **0.3400**)

Ranking of tested nearby variations:
1. `candidate_variant_a` — median peripyramidal reactive-astrocyte area: partial r **-0.3755**, selection **0.3755**, LOO **0.0995**
2. `candidate_variant_b` — 75th percentile peripyramidal reactive-astrocyte area: partial r **-0.2306**, selection **0.2306**, LOO **0.0342**

So the robust median beat the upper-tail summary, but both variants showed weak out-of-sample behavior.

## Findings
1. **What worked and why**
   - Restricting to **CA1 Reactive Astrocytes** in the **CA1 pyramidal-neuron niche** preserved the same biologically focused lead as earlier rounds while changing the donor scalar from burden/crowding to **hypertrophy severity**.
   - The best donor scalar was the **median contour area** of reactive astrocytes whose centroids lay within **30 µm** of any CA1 pyramidal neuron. This likely works better than the 75th percentile because it captures a donor-wide shift toward larger perineuronal reactive astrocytes rather than a few giant cells.
   - Coverage stayed high because every donor had qualifying cells, with a median of roughly a few hundred peripyramidal reactive astrocytes per slide.

2. **What failed and why**
   - The signal did **not** generalize well: LOO predictive r was only **0.0995** despite an in-sample partial r of **-0.3755**.
   - The upper-tail alternative (`candidate_variant_b`) underperformed, suggesting that the most extreme hypertrophic cells are too sparse or donor-specific to be a stable scalar on this cohort.
   - This family looks partly redundant with the accepted CA1 reactive-astrocyte panel axis: it tracks the same peripyramidal astrocyte process, but size alone does not yet produce a robust new predictive signal beyond the earlier burden/exposure/crowding summaries.

3. **Error pattern: which donors are consistently wrong and what they share**
   - The most influential donor was **H20.33.018**; leaving it out shifted the partial correlation by about **0.3400**, making it the only unstable donor under the round threshold.
   - The largest underpredictions included **H20.33.018**, **H20.33.038**, and **H20.33.046**. These donors share relatively severe pathology labels in the available confounds (mostly **Braak 5–6**, **CERAD 3**) and more negative memory decline than the feature alone predicted.
   - The largest overpredictions included **H21.33.040**, **H21.33.035**, and **H20.33.024**: donors with moderate-to-high pathology and elevated peripyramidal astrocyte size but less severe observed decline than expected.
   - Interpretable pattern: peripyramidal reactive-astrocyte hypertrophy appears to track a disease-related niche state, but once pathology is already high, **size alone saturates** and does not separate the most impaired from the merely pathology-heavy donors.

## Rationale
The winning approach is biologically coherent because it stays anchored to a clear cell population and niche: **reactive astrocytes around CA1 pyramidal neurons**, where astrocytic hypertrophy could reflect localized glial activation around vulnerable excitatory neurons. It beat the nearby 75th-percentile alternative because the median is more robust to a few extreme contours and better represents a donor-level shift in the typical perineuronal reactive-astrocyte state.

That said, the weak LOO behavior suggests this feature is more likely a **refinement of the existing CA1 reactive-astrocyte axis** than a clearly additive new panel member. It is biologically sensible, but not yet compelling as a replacement or addition on predictive grounds alone.

## Interpretation
Biologically, the signal seems to mean: donors with worse memory decline tend to show **larger, hypertrophic reactive astrocytes in the CA1 pyramidal-neuron neighborhood**.

- **Population:** CA1 Reactive Astrocytes
- **Niche:** within 30 µm of a CA1 Pyramidal Neuron centroid
- **Feature summary:** donor-level median contour area in µm²
- **Simplest observable pattern:** a thicker, enlarged ring of reactive astrocyte profiles clustered around the CA1 pyramidal layer

This is a morphology-based readout of local astrocyte activation severity rather than pure burden.

## Next
Stay in the same niche and test a **hypertrophy prevalence** sweep rather than another tail statistic: for example, the fraction of CA1 peripyramidal reactive astrocytes above a fixed large-area threshold or above a donor-external reference threshold. That directly addresses the error pattern seen here, where the median beat the 75th percentile, implying donor-wide hypertrophy burden may be more stable than rare giant-cell extremes.
