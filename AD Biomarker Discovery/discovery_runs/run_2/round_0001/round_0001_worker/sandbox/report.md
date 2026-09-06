## Summary
Tested the CA1 pyramidal neuron burden family with two local normalization variants; `ca1_pyramidal_fraction` won with selection score 0.2940.

## Metrics
Winning variation: `ca1_pyramidal_fraction`.
- IS partial r: 0.2940
- Selection score: 0.2940
- LOO predictive r: 0.0976
- IS-LOO gap: 0.1964
- Penalty: 0.0964
- Adjusted score: 0.1976
- Analyzable donors: 35/35

Ranking of tested variations:
1. `ca1_pyramidal_fraction`: partial r 0.2940, selection score 0.2940, LOO predictive r 0.0976.
2. `ca1_pyramidal_density`: partial r 0.2122, selection score 0.2122, LOO predictive r 0.0246.

## Findings
1. What worked and why: The fraction variant worked best because it normalizes CA1 pyramidal burden against the total classified CA1 cellular compartment, reducing slide-to-slide differences in annotation size and overall tissue amount.
2. What failed and why: The density variant underperformed because absolute counts per annotated area appear more sensitive to regional outline size, section thickness, and local packing differences that are not specific to selective CA1 neuron loss.
3. Error pattern: Worst absolute LOO errors were H20.33.018 (|err|=0.191), H21.33.040 (|err|=0.137), H21.33.023 (|err|=0.134), H21.33.004 (|err|=0.101), H21.33.046 (|err|=0.098); among these donors, sex `Male` (4/5).

## Rationale
A compositional loss signal is biologically coherent in CA1: as pyramidal neurons are depleted or replaced by glial and other non-pyramidal cells, the pyramidal share of CA1 should fall. This beat raw density because the nearby alternative retained extra geometric variance from the annotation footprint.
Because this is round 1 with no accepted panel, the winner is a plausible seed biomarker rather than evidence of additivity beyond an existing panel.

## Interpretation
Biologically, the signal looks like CA1 pyramidal neuron preservation versus depletion within hippocampal memory circuitry.
- Population: `Pyramidal Neuron`
- Niche/region: `CA1`
- Donor-level scalar: share of all CA1 classified cells labeled Pyramidal Neuron
- Simplest observable tissue pattern: a thinner CA1 pyramidal band relative to the total local CA1 cell population

## Next
Next, keep the CA1 pyramidal-burden direction but test denominators that are closer to the CA1 circuit itself, for example CA1 pyramidal cells divided by all CA1 neuronal cells rather than all classified CA1 cells.
