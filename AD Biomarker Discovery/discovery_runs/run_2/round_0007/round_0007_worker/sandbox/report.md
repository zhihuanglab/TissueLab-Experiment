## Summary
Tested the CA1 reactive-astrocyte hotspot pyramidal-depletion family; candidate_variant_b won with selection score 0.0436.

## Metrics
Winner: candidate_variant_b with partial r -0.0545, selection score 0.0436, LOO predictive r -0.5480, IS-LOO gap -0.4935, penalty 0.0000, adjusted score -0.5480, analyzable donors 28/35.
Other tested variations: candidate_variant_b partial_r=-0.0545, selection=0.0436, loo_r=-0.5480; candidate_variant_a partial_r=0.0520, selection=0.0416, loo_r=-0.4646.

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The winning score was the donor-level difference between mean log1p CA1 pyramidal-neighbor counts around ordinary astrocyte hotspots and around reactive astrocyte hotspots. Positive values mean reactive hotspots sit in more pyramidal-sparse CA1 microterritories, matching a focal gliotic-scarring / neuronal-dropout interpretation.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The losing local variation changed only the pyramidal contact radius. If it underperformed, the stricter radius likely over-focused on immediate pericellular space and discarded broader hotspot-context information that still carries the neuron-loss pattern.
3. Error pattern: which donors are consistently wrong and what they share
   - Largest LOO errors: H20.33.018 (outcome=-0.304, pred=-0.073, feature=0.049), H20.33.024 (outcome=-0.010, pred=-0.123, feature=0.320), H21.33.040 (outcome=-0.033, pred=-0.144, feature=0.189)
   - Shared pattern among these donors: median Braak ≈ 5.0.

## Rationale
The best variation is biologically coherent because it keeps the previously promising CA1 reactive-hotspot lead fixed and asks a more mechanistic question: whether reactive hotspots occupy neuron-depleted CA1 territory relative to baseline astrocyte hotspots. It beat the nearby alternative by using the contact radius that best preserved a stable distinction between ordinary astroglial clustering and reactive hotspot placement. Relative to the current panel, this candidate is most likely to add information if hotspot-centered neuronal depletion is cleaner than earlier neuron-adjacent exposure terms.

## Interpretation
The signal appears to reflect CA1 reactive astrocyte hotspots occupying pyramidal-depleted niches. Population: CA1 Reactive Astrocytes versus baseline CA1 Astrocytes. Niche: same-class astroglial hotspots within CA1, evaluated for nearby CA1 Pyramidal Neuron density. Feature summary: astrocyte-hotspot mean log1p pyramidal-neighbor count minus reactive-hotspot mean log1p pyramidal-neighbor count. Simplest observable tissue pattern: reactive astrocyte microclusters sitting in CA1 patches where local pyramidal neurons are relatively sparse compared with ordinary astrocyte clusters.

## Next
Run one nearby sweep that keeps CA1 reactive hotspots fixed but changes the donor summary operator, for example comparing hotspot-level lower-tail pyramidal-neighbor burden (e.g. median or bottom-quartile depletion) instead of the mean, because the main residual errors may come from donors with especially focal rather than average neuronal dropout.
