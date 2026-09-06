## Summary
Tested the CA1 reactive-astrocyte cuff severity around CA1 pyramidal neurons within 35 µm family, and candidate_variant_a won with selection score 0.5044.

## Metrics
- Best variation: `candidate_variant_a` (`ca1_pyramidal_cuffed_fraction_r35um_ge2`)
- IS partial r: `-0.5044`
- Selection score: `0.5044`
- LOO predictive r: `0.3412`
- IS-LOO Gap: `0.1633` `(penalty=0.0133)`
- Adjusted Score: `0.3279`
- Other tested variations: candidate_variant_b selection=0.4997, partial_r=-0.4997, loo_r=0.3540

## Findings
1. What worked and why (tie to the biological meaning of the target)  
   The winning feature was `ca1_pyramidal_cuffed_fraction_r35um_ge2`, the donor-level fraction of CA1 pyramidal neurons that sit within a 35 µm cuff of at least 2 CA1 reactive astrocytes. This worked best because it keeps the niche that already looked promising in the accepted panel—reactive astrocytes wrapped around CA1 pyramidal neurons—but asks for more than a single nearby astrocyte, which should better capture true perineuronal reactive cuffs rather than incidental proximity.
2. What failed and why (specific to the chosen hypothesis and what went wrong)  
   The nearby stricter threshold lost to the winner: candidate_variant_b selection=0.4997, partial_r=-0.4997, loo_r=0.3540. Requiring too many nearby reactive astrocytes appears to make the cuffing event too sparse, so the feature starts to emphasize only the most extreme local lesions and loses donor-level dynamic range.
3. Error pattern: which donors are consistently wrong and what they share  
   The largest absolute LOO errors were H20.33.018 (Braak 6, Female), H21.33.018 (Braak 2, Female), H21.33.023 (Braak 4, Male). The biggest negative errors (observed memory slope lower than predicted) were H21.33.043 (Braak 2, Female), H21.33.006 (Braak 4, Male), H21.33.007 (Braak 5, Female), while the biggest positive errors were H20.33.018 (Braak 6, Female), H21.33.018 (Braak 2, Female), H21.33.023 (Braak 4, Male). These misfit donors suggest the cuffing signal is strongest when CA1 neuron-centered reactive astrocytosis is a major driver, but some donors likely carry additional decline mechanisms not fully captured by this single niche summary.

## Rationale
The best variation is biologically coherent because it encodes a specific microenvironment: CA1 pyramidal neurons surrounded by multiple reactive astrocytes inside a short perineuronal radius. That is closer to an interpretable injury-response cuff than a whole-region astrocyte burden. It beat the nearby alternative because the winning threshold preserved enough prevalence across donors to remain measurable, while the stricter threshold likely over-focused on rare extremes. Relative to the current panel, it is plausibly additive only if the exact multi-astrocyte severity threshold captures a sharper state than the existing any-cuffing feature; otherwise it may mostly act as a refinement or possible replacement rather than a fully new axis.

## Interpretation
The signal appears to mean the burden of CA1 pyramidal neurons sitting inside compact reactive-astrocyte cuffs.  
Population: CA1 Pyramidal Neuron targets with CA1 Reactive Astrocyte neighbors.  
Niche: a 35 µm perineuronal cuff centered on pyramidal-neuron centroids in CA1.  
Feature summary: donor-level fraction of CA1 pyramidal neurons with at least 2 reactive astrocytes inside that radius.  
Simplest observable pattern: more CA1 pyramidal neurons visibly ringed by several reactive astrocytes corresponds to worse memory-related slope.

## Next
Run one more local sweep that keeps the same CA1 pyramidal-centered cuff niche and 35 µm radius but tests a mild intensity refinement such as reactive-neighbor count percentiles or mean count among already-cuffed neurons, because this round tells us thresholding matters and the main loss mode is likely oversparsifying the severe-end definition.
