## Summary
Tested the ca1_peripyramidal_reactive_astro_coverage family in CA1, and coverage_60px won with selection_score 0.5176.

## Metrics
Winning variation: coverage_60px (60 px).
- IS partial r: -0.5176
- Selection score: 0.5176
- LOO predictive r: 0.3705
- IS-LOO Gap: 0.1471 (penalty=0.0000)
- Adjusted Score: 0.5176
- Analyzable donors: 35/35

Other tested variations:
- coverage_80px: partial_r=-0.5156, selection_score=0.5156, loo_predictive_r=0.3700, adjusted_score=0.5156

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The winner works best when the signal is expressed as a neuron-anchored CA1 exposure metric: among surviving CA1 pyramidal neurons, what fraction sits within 60 px of at least one reactive astrocyte. This sharpens the biology from bulk cell composition toward a local neuron-glia interaction state that plausibly tracks stressed or remodeled CA1 tissue relevant to memory decline.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The broader nearby radius was slightly weaker, which suggests that expanding the neighborhood starts to dilute the specifically perineuronal reactive-astrocyte signal with more background CA1 astrocyte burden.
3. Error pattern: which donors are consistently wrong and what they share
- H21.33.018: outcome=-0.1440, predicted=0.0050, ca1_pyramidal_reactive_astro_coverage_60px=0.0634, CA1 pyramidal=1783, CA1 reactive astrocyte=345
- H20.33.018: outcome=-0.3042, predicted=-0.1617, ca1_pyramidal_reactive_astro_coverage_60px=0.3912, CA1 pyramidal=1020, CA1 reactive astrocyte=1161
- H21.33.023: outcome=-0.1546, predicted=-0.0284, ca1_pyramidal_reactive_astro_coverage_60px=0.0856, CA1 pyramidal=1881, CA1 reactive astrocyte=449
- H20.33.046: outcome=-0.2017, predicted=-0.0947, ca1_pyramidal_reactive_astro_coverage_60px=0.1156, CA1 pyramidal=2058, CA1 reactive astrocyte=482
- H21.33.006: outcome=-0.0540, predicted=-0.1553, ca1_pyramidal_reactive_astro_coverage_60px=0.2482, CA1 pyramidal=1378, CA1 reactive astrocyte=1010

## Rationale
The best variation is biologically coherent because it conditions first on the vulnerable CA1 pyramidal neuron population and then measures local reactive astrocyte proximity around those neurons, rather than counting reactive astrocytes globally. That makes the scalar interpretable as peripyramidal reactive astrocyte coverage. It beat the nearby alternative because the winning radius stayed tighter around the perineuronal zone and avoided diluting the signal with more diffuse background CA1 reactive astrocyte burden. Relative to the accepted panel, this candidate is plausibly somewhat additive if it captures conditional local exposure beyond overall CA1 pyramidal loss, reactive astrocyte burden, and the broader mixed-cell niche fraction already in panel review.

## Interpretation
Biologically, the signal appears to represent how strongly surviving CA1 pyramidal neurons are embedded in a reactive astrocyte-rich local niche. The population is CA1 Pyramidal Neuron, the niche partner is CA1 Reactive Astrocyte, the donor-level scalar is the fraction of CA1 pyramidal neurons with at least one reactive astrocyte within 60 px, and the simplest observable tissue pattern is a higher share of CA1 pyramidal neurons with nearby reactive astrocyte neighbors crowding the local perineuronal space.

## Next
Run the next local sweep around the same neuron-anchored CA1 reactive-astrocyte exposure idea but vary only the donor-level aggregation, for example comparing binary coverage to mean nearest-reactive-astrocyte distance or a soft distance-weighted coverage score, because the error pattern suggests the family has signal but may be losing information by thresholding all covered neurons equally.
