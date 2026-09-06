## Summary
Tested the CA1 pyramidal-near-reactive-astrocyte niche family with a 30 um vs 50 um radius sweep; pyramidal_near_reactiveastro_30um won with selection_score 0.5191.

## Metrics
- Winning variation: pyramidal_near_reactiveastro_30um
- IS partial r: -0.5191
- Selection score: 0.5191
- LOO predictive r: 0.3727
- IS-LOO Gap: 0.1464 (penalty=0.0000)
- Adjusted Score: 0.3727
- Other tested variations: pyramidal_near_reactiveastro_50um partial_r=-0.4987, selection_score=0.4987, loo_r=0.3534

## Findings
1. What worked and why: The strongest signal came from the 30 um peri-pyramidal niche in CA1, where the donor scalar is the fraction of CA1 pyramidal neurons that have at least one nearby reactive astrocyte. That beats the broader local alternative, which suggests the biomarker is tied to a fairly immediate neuron-adjacent reactive astrocyte pattern rather than diffuse CA1 gliosis alone.
2. What failed and why: The broader 50 um neighborhood underperformed because it likely dilutes the perineuronal niche with more background reactive astrocytes, making the feature look more like general CA1 reactive astrocyte burden and less like a sharply neuron-adjacent injury pattern.
3. Error pattern: Largest absolute LOO errors were for H21.33.018 (outcome=-0.144, pred=0.005, 30um=0.063, 50um=0.164, CA1 pyramidal=1783, CA1 reactive astro=345), H20.33.018 (outcome=-0.304, pred=-0.162, 30um=0.391, 50um=0.715, CA1 pyramidal=1020, CA1 reactive astro=1161), H21.33.023 (outcome=-0.155, pred=-0.029, 30um=0.086, 50um=0.230, CA1 pyramidal=1881, CA1 reactive astro=449).

## Rationale
This winner is biologically coherent because it keeps the round-1 CA1 reactive astrocyte signal but sharpens it to the specific niche around CA1 pyramidal neurons. The summary is simple and auditable: among all CA1 pyramidal neurons, what fraction lies within 30 um of a reactive astrocyte? It beat the nearby radius alternative because that distance best balanced specificity to a neuron-adjacent niche against donor-level stability. Because it uses the same CA1 reactive astrocyte biology as the accepted panel member, it is more likely to refine or replace that member than to add wholly orthogonal information.

## Interpretation
The signal appears to represent CA1 reactive astrocyte engagement in the immediate peri-pyramidal niche. Population: CA1 pyramidal neurons and CA1 reactive astrocytes. Niche: reactive astrocytes within 30 um of pyramidal somata. Feature summary: donor-level fraction of CA1 pyramidal neurons that are niche-positive. Simplest observable tissue pattern: CA1 fields where many pyramidal neurons sit directly alongside reactive astrocytes track worse memory decline.

## Next
If this niche winner is kept, the next local sweep should stay in CA1 pyramidal–reactive astro biology but test whether normalizing by overall CA1 reactive astro burden or restricting to a very tight inner band improves additivity over the current panel.
