## Summary
One sentence: tested the CA1 pyramidal reactive-astrocyte cuffing fraction family, r35_ge1 won, and its local winner selection score was 0.5224.

## Metrics
Winning variation: r35_ge1.
- IS partial r: -0.5224
- Selection score: 0.5224
- LOO predictive r: 0.3768
- IS-LOO Gap: 0.1456 (penalty=0.1456)
- Adjusted Score: 0.3768

Other tested variations ranked by local single-feature signal:
- r35_ge1 (selection=0.5224, partial_r=-0.5224, loo_r=0.3768); r35_ge2 (selection=0.5044, partial_r=-0.5044, loo_r=0.3412); r35_ge3 (selection=0.4997, partial_r=-0.4997, loo_r=0.3540)

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The neuron-centric cuffing summary worked best when defined as r35_ge1, meaning the donor-level signal is strongest when asking how broadly CA1 pyramidal neurons are exposed to nearby reactive astrocytes, rather than only whether a minimal contact exists. This keeps the population fixed to CA1 pyramidal neurons and sharpens the niche to a 35 µm perineuronal neighborhood.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The nearby thresholds lost signal because they were either too permissive (r35_ge2) and likely collapsed donors toward saturation, or too strict (r35_ge3) and likely discarded informative moderate cuffing events. The family appears sensitive to how much reactive-astrocyte occupancy is required before a neuron counts as exposed.
3. Error pattern: which donors are consistently wrong and what they share
   - Largest LOO errors were:
- H21.33.018: outcome=-0.144, predicted=0.010, ca1_pyramidal_cuffed_fraction_r35um_ge1=0.081
- H20.33.018: outcome=-0.304, predicted=-0.162, ca1_pyramidal_cuffed_fraction_r35um_ge1=0.486
- H21.33.023: outcome=-0.155, predicted=-0.025, ca1_pyramidal_cuffed_fraction_r35um_ge1=0.113
- H20.33.046: outcome=-0.202, predicted=-0.098, ca1_pyramidal_cuffed_fraction_r35um_ge1=0.168
- H21.33.006: outcome=-0.054, predicted=-0.151, ca1_pyramidal_cuffed_fraction_r35um_ge1=0.307

## Rationale
Why the best variation's approach is biologically coherent and why it beat the nearby alternatives.
- CA1 pyramidal neurons are the hippocampal population most directly tied to memory circuitry, so measuring how often they sit inside a local ring of reactive astrocytes is a biologically coherent niche summary.
- This beat the nearby alternatives because the winning threshold best balanced prevalence and specificity: enough reactive astrocytes to indicate a real cuffing state, but not so many that only extreme cases contribute.
- Relative to the current panel, it is plausibly additive because it converts an astrocyte-centered reactivity niche into a neuron-exposure niche.

## Interpretation
State this when it materially helps explain the biomarker.
- The signal appears to mean that donors with worse memory decline have a larger fraction of CA1 pyramidal neurons locally surrounded by reactive astrocytes.
- Population: CA1 Pyramidal Neuron, with neighboring CA1 Reactive Astrocyte.
- Niche: a 35 µm centroid-defined perineuronal neighborhood inside CA1.
- Feature summary: donor-level fraction of CA1 pyramidal neurons whose local reactive-astrocyte neighbor count crosses the winning threshold.
- Simplest observable pattern: broader reactive astrocyte cuffing around CA1 pyramidal neurons.

## Next
One specific suggestion for the next local sweep based on the error pattern and which nearby variations won or lost.
- Keep the same CA1 pyramidal / reactive astrocyte niche but sweep radius around the winning threshold setting (for example 25, 45, and 55 µm with the winning minimum-neighbor rule fixed) to test whether the current errors reflect an overly narrow spatial scale rather than the exposure definition itself.
