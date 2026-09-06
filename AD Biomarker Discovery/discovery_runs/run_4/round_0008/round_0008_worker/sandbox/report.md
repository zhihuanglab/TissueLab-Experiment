## Summary
Tested the `ca1_reactive_oligodendrocyte_gap` family in CA1 reactive-astrocyte-conditioned pyramidal niches; `candidate_variant_b` (80 px) won with selection_score=0.1351.

## Metrics
Best variation: `candidate_variant_b` using an 80 px radius.
- IS partial r: 0.1351
- Selection score: 0.1351
- LOO predictive r: -0.1228
- IS-LOO Gap: 0.2579 (penalty=0.0539)
- Adjusted Score: -0.1767
- n_analyzable / n_total: 35 / 35

Other tested variation:
- `candidate_variant_a` (60 px): partial_r=0.0922, selection_score=0.0922, loo_predictive_r=-0.2103, adjusted_score=-1.0000

So the broader local niche won inside the family, but both variants were weak out of sample and the family does not look panel-ready in its current form.

## Findings
1. What worked and why: The best local signal came from CA1 pyramidal neurons already exposed to reactive astrocytes, then asking whether oligodendrocytes were absent nearby. That is biologically coherent because it isolates a stressed neuronal niche where injury-reactive glia are present but local myelin/support-cell coverage is sparse.
2. What failed and why: The whole family was unstable out of sample. The 60 px version was too narrow and likely missed the broader support-cell desert around reactive astrocyte niches; the 80 px version improved the in-sample partial correlation slightly, but the negative LOO r shows that local oligodendrocyte absence is not yet robust enough as a standalone donor scalar.
3. Error pattern: The biggest misses were H20.33.018, H21.33.023, H21.33.040, H20.33.024, and H20.33.037. These are mostly older donors (79-100 years) with Braak 4-6 pathology. The model underpredicted several fast-decline donors with low gap fractions (for example H20.33.018 and H20.33.037), while it overpredicted decline in milder-outcome donors with only moderate gap fractions (for example H21.33.040 and H20.33.024). That pattern suggests the feature tracks one aspect of CA1 niche injury but misses additional severity axes already mixed into the accepted panel.

## Rationale
The winning variation is biologically coherent because it conditions on the reactive astrocyte-exposed CA1 pyramidal population instead of measuring oligodendrocytes globally. That denominator asks a specific mechanistic question: when CA1 neurons sit inside a reactive glial neighborhood, how often is oligodendrocyte support locally absent?

It beat the 60 px alternative because oligodendrocyte depletion appears to operate on a slightly broader spatial scale than the tightest reactive contact zone. Even so, the family seems partly redundant with the current CA1 injury panel: once CA1 neuron loss, reactive enrichment, and reactive-pyramidal niche injury are already represented, a binary oligodendrocyte-gap fraction does not add a stable new axis by itself.

## Interpretation
Biologically, the signal seems to mean that worse memory decline is associated with CA1 neuronal neighborhoods where reactive astrocytes are close to pyramidal neurons but oligodendrocytes are locally sparse.

- Population: CA1 Pyramidal Neuron
- Niche: CA1 reactive-astrocyte-exposed pyramidal neighborhoods
- Feature summary: fraction of reactive-exposed CA1 pyramidal neurons with zero oligodendrocyte neighbors within 80 px
- Simplest observable pattern: reactive astrocyte-wrapped CA1 pyramidal zones with few nearby oligodendrocytes / weak local myelin-support presence

## Next
Stay in the same CA1 pyramidal/reactive/oligodendrocyte family, but replace the binary zero-neighbor fraction with a graded local-support score at the winning scale (around 80 px), such as mean oligodendrocyte neighbor count or distance-to-nearest-oligodendrocyte among reactive-exposed pyramidal neurons.
