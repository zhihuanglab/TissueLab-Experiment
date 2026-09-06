## Summary
One sentence: tested the CA1 peripyramidal reactive-astrocyte hypertrophy within 35 µm of CA1 pyramidal neurons family, candidate_variant_a won (Median contour area of CA1 reactive astrocytes within 35 µm of any CA1 pyramidal neuron), and the local winner selection score was 0.3165.

## Metrics
Winner: candidate_variant_a
- IS partial r: -0.3165
- Selection score: 0.3165
- LOO predictive r: -0.0289
- IS-LOO Gap: 0.2876 (penalty=0.1376)
- Adjusted Score: -0.1087

Other tested variations: candidate_variant_b (selection=0.2523, partial_r=-0.2523, loo=0.0853)

## Findings
1. What worked and why: The winning summary favored the central tendency of qualifying reactive-astrocyte size, suggesting that broad niche-wide hypertrophy is more reproducible than focusing only on the most extreme cells.
2. What failed and why: The losing variation likely over- or under-emphasized the hypertrophic tail. Within this fixed niche, changing only the donor-level area summary moved signal less than earlier count-based localization steps, suggesting morphology adds only a modest refinement on top of the already-established cuffing niche.
3. Error pattern: the largest LOO errors were H20.33.018 (|e|=0.215), H20.33.038 (|e|=0.122), H21.33.035 (|e|=0.113), H21.33.040 (|e|=0.112), H20.33.046 (|e|=0.107); 3/5 of the largest errors are on donors with more negative-than-median memory slope.

## Rationale
The best variation is biologically coherent because it stays inside the same CA1 peripyramidal reactive-astrocyte niche that drove rounds 3-5, but changes the donor summary from abundance to morphology. Reactive astrocyte contour area is a plausible proxy for hypertrophic activation state on LFB tissue, and the winner beat the nearby alternative because its chosen summary better balanced niche specificity against donor-level robustness. Because the current accepted panel already contains broad reactivity and CA1 cuffing features, this morphology summary is biologically coherent but may be partially redundant unless its hypertrophy readout captures severity beyond simple neighbor counts.

## Interpretation
The signal appears to reflect a donor-level shift toward generally enlarged reactive astrocytes in the CA1 peripyramidal niche. Population: CA1 Reactive Astrocytes. Niche: within 35 µm of CA1 Pyramidal Neurons. Feature summary: median qualifying reactive-astrocyte area. Simplest observable pattern: larger swollen reactive astrocyte profiles clustered directly around CA1 pyramidal neurons track worse memory decline.

## Next
Next local sweep: keep the same CA1 peripyramidal reactive-astrocyte niche but test whether hypertrophy is most informative when conditioned on reactive-cuff burden (for example, median or p75 area only among donors' cuffed pyramidal neighborhoods with >=2 nearby reactive astrocytes), since this round asks morphology alone to add on top of already-strong count-based winners.
