## Summary
One sentence: tested the CA1 peripyramidal astrocytic reactivity family, peripyramidal_reactivity_r35um won, and its local winner selection score was 0.4941.

## Metrics
Winning variation `peripyramidal_reactivity_r35um` used feature column `ca1_peripyramidal_reactivity_r35um` with partial r -0.4941, selection score 0.4941, LOO predictive r 0.3441, IS-LOO gap 0.1500, penalty 0.1500, and adjusted score 0.3441.
- `peripyramidal_reactivity_r25um`: partial r -0.4759, selection score 0.4759, LOO predictive r 0.3208.

## Findings
1. What worked and why (tie to the biological meaning of the target)  
   The winner sharpened the accepted astrocytic-reactivity biology to a CA1 neuron-adjacent niche: reactive astrocytes among CA1 astrocyte-lineage cells lying within 35 µm of a CA1 pyramidal neuron. This directly targets peripyramidal gliosis around the neuron population most plausibly linked to memory decline.
2. What failed and why (specific to the chosen hypothesis and what went wrong)  
   The nearby radius alternative underperformed, suggesting the signal is not improved by simply broadening the neighborhood. If the broader radius lost, the niche is probably fairly tight; if the tighter radius lost, the signal likely needs a slightly more permissive peripyramidal field.
3. Error pattern: which donors are consistently wrong and what they share  
The largest LOO errors are concentrated in Male donors with median Braak 4.0 and a median peripyramidal denominator of 512; these appear to be donors where CA1 reactive fraction alone does not fully track outcome severity.
- H20.33.018: outcome -0.304, predicted -0.141, |error| 0.163, sex Female, Braak 6, CERAD 3, denominator 669.
- H21.33.023: outcome -0.155, predicted 0.006, |error| 0.161, sex Male, Braak 4, CERAD 0, denominator 512.
- H21.33.006: outcome -0.054, predicted -0.190, |error| 0.136, sex Male, Braak 4, CERAD 1, denominator 496.

## Rationale
The best variation is biologically coherent because it asks whether the CA1 astrocyte-reactivity signal is concentrated specifically around CA1 pyramidal neurons rather than diffusely across all of CA1. That beat the nearby alternative because the winning radius better matched the scale at which astrocyte activation appears to covary with donor memory-decline burden. Since the current accepted panel already contains a region-wide astrocytic reactivity member, this niche-refined candidate is likely partly redundant but could still add value if the evaluator finds that neuron-adjacent gliosis carries cleaner or slightly more specific signal.

## Interpretation
The signal seems to represent CA1 peripyramidal gliosis: the population is CA1 astrocyte-lineage cells, the niche is immediate proximity to CA1 pyramidal neurons, the donor-level summary is the reactive fraction among those niche astrocytes, and the simplest observable pattern is more reactive astrocytes crowding the CA1 pyramidal field in donors with worse memory decline.

## Next
Run one more local sweep around the winning niche scale by testing a slightly tighter or slightly broader CA1 pyramidal-neighborhood threshold and, if needed, a denominator-stabilized variant such as a shrinkage-adjusted reactive fraction.
