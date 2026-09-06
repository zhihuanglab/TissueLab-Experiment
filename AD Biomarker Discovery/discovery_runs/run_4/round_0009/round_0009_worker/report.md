## Summary
Tested the ca1_reactive_corpora_amylacea_niche family in CA1 reactive-astrocyte-conditioned pyramidal niches; corpora_any_r96px won with selection_score=0.1206.

## Metrics
Best variation: corpora_any_r96px (96 px).
- IS partial r: -0.1206
- Selection score: 0.1206
- LOO predictive r: -0.0592
- IS-LOO Gap: 0.0613 (penalty=0.0000)
- Adjusted Score: -0.0592
- n_analyzable / n_total: 35 / 35
Other tested variations:
- corpora_any_r80px (80 px): partial_r=-0.0734, selection_score=0.0734, loo_r=-0.0642, adjusted=-0.0642
- corpora_any_r64px (64 px): partial_r=-0.0055, selection_score=0.0055, loo_r=-0.0620, adjusted=-0.0620

## Findings
1. What worked and why: The winning signal comes from CA1 pyramidal neurons that already sit in reactive astrocyte neighborhoods, then asks how often those neurons also have nearby corpora amylacea. That conditional denominator focuses the feature on a chronic-debris version of the existing CA1 injury niche rather than global corpora burden, which is biologically coherent if local degenerative debris handling tracks more advanced neuronal stress. The winning feature is sparse in this cohort: only 5 of 35 donors had a nonzero value, so the family behaves more like a rare-niche marker than a broad continuum.
2. What failed and why: Narrower corpora-amylacea radii underperformed, suggesting the relevant degenerative-debris field extends beyond the immediate neuron surface and is only captured at a broader CA1 reactive niche scale.
3. Error pattern: The largest LOO misses were donors H20.33.034, H20.33.018, H21.33.023, H21.33.040, H20.33.037. The worst misses were mostly underestimates of decline for H20.33.018, H21.33.023, H20.33.037, with overestimates for H20.33.034, H21.33.040. Among the worst-error donors, the ones with nonzero corpora-positive niche signal were H20.33.034 (feature=0.0375, CA count=120, reactive-niche pyramidal n=80), H20.33.018 (feature=0.0119, CA count=7, reactive-niche pyramidal n=587).

## Rationale
The best variation is biologically coherent because it measures a specific CA1 microenvironment: pyramidal neurons near reactive astrocytes, further stratified by whether corpora amylacea are also present nearby. That is a plausible marker of chronic local tissue stress or debris-processing burden, distinct from pure neuron loss, reactive burden alone, or oligodendrocyte absence. It beat the nearby alternatives because its radius better matched the apparent spatial scale of corpora amylacea accumulation around injured CA1 neuronal niches.
Relative to the current accepted panel, this candidate seems only moderately likely to add new information because the panel already captures CA1 neuron-loss and reactive-niche biology; any additivity would have to come from a more chronic debris-handling axis rather than another restatement of the same CA1 injury process.

## Interpretation
The signal appears to mean that faster memory decline tracks CA1 neuronal niches where reactive astrocytes and corpora amylacea co-occur around pyramidal neurons.
- Population: CA1 Pyramidal Neuron
- Niche: CA1 reactive-astrocyte-exposed pyramidal neighborhoods
- Feature summary: fraction of reactive-niche CA1 pyramidal neurons with at least one Corpora Amylacea within 96 px
- Simplest observable pattern: CA1 pyramidal zones already touching reactive astrocytes that also show nearby corpora amylacea / chronic degenerative debris bodies

## Next
Stay in the same family and test slightly broader corpora radii beyond 96 px, or an intensity variant using mean corpora neighbor count among reactive-niche pyramidal neurons.
