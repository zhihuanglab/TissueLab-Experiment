## Summary
Tested the EC reactive-astrocyte-proximal pyramidal lower-tail log-area family in EC across 80 px and 60 px reactive-astrocyte proximity radii; `candidate_variant_a` won with selection_score `0.2205`.

## Metrics
- Best variation: `candidate_variant_a` (`ec_reactive_astro_proximal_pyramidal_lower_tail_area_80px`)
- IS partial r: `0.2269`
- Selection score: `0.2205`
- LOO predictive r: `0.0326`
- IS-LOO Gap: `0.1943` `(penalty=0.0221)`
- Adjusted Score: `0.0105`
- Coverage: `34/35` donors
- Other tested variations: `candidate_variant_b` selection_score `0.1526`, partial_r `0.1571`, LOO `-0.1139`

## Findings
1. What worked and why (tie to the biological meaning of the target)  
   The best signal came from EC pyramidal neurons that sit within `80` px of EC reactive astrocytes, summarized by the 10th percentile of `log1p(contour area)`. That focuses the biomarker on the lower tail of neuronal size inside a glial-reactive niche rather than on global EC composition, which is biologically coherent for memory decline because it reads out a locally stressed excitatory-neuron subpopulation in entorhinal cortex.
2. What failed and why (specific to the chosen hypothesis and what went wrong)  
   `candidate_variant_b` (radius 60 px) reached selection_score `0.1526` with partial_r `0.1571` and LOO predictive r `-0.1139`. The losing radius likely either diluted the pericellular niche by admitting more distant pyramidal neurons or over-tightened the niche and reduced stability, depending on whether it was broader or tighter than the winner.
3. Error pattern: which donors are consistently wrong and what they share  
   Largest underpredictions were H20.33.018 (err=-0.2020, braak=6, cerad=3, sex=Female), H21.33.023 (err=-0.1453, braak=4, cerad=0, sex=Male), H20.33.038 (err=-0.1074, braak=5, cerad=3, sex=Female). Largest overpredictions were H21.33.040 (err=+0.1202, braak=5, cerad=3, sex=Male), H20.33.024 (err=+0.0785, braak=5, cerad=1, sex=Male), H21.33.035 (err=+0.0741, braak=5, cerad=2, sex=Female). Mean Braak values were `5.00` for underpredictions and `5.00` for overpredictions, while mean proximal pyramidal counts were `1110.3` and `386.0` respectively. That suggests remaining errors are not just sparse-count failures but donor-specific mismatch between EC niche morphology and memory decline.

## Rationale
The winning variation is biologically coherent because it isolates one population (EC pyramidal neurons), one niche sharpener (nearest EC reactive astrocytes within a short Euclidean radius), and one donor scalar (the lower tail of `log1p(area)`). That combination should enrich for small, morphologically contracted pyramidal neurons in a locally reactive glial microenvironment, a plausible tissue correlate of neurodegenerative stress relevant to memory. It beat the nearby alternative because `80` px appears to strike the better balance between niche specificity and donor-level stability on this cohort. Relative to the accepted CA1-focused panel, this EC transfer is mechanistically aligned but regionally distinct, so it may still add information if the evaluator finds it non-redundant.

## Interpretation
The signal seems to mean that donors with worse memory decline have a lower lower-tail area among EC pyramidal neurons that lie close to reactive astrocytes. Population: `Pyramidal Neuron`. Niche: EC cells within `80` px of an EC `Reactive Astrocyte`. Feature summary: 10th percentile of `log1p(contour area)`. Simplest observable pattern: in entorhinal cortex, the reactive-astrocyte-adjacent pyramidal neurons look smaller/atrophic in worse-decline donors.

## Next
Keep EC fixed but sweep the lower-tail summary itself next: compare the 5th, 10th, and 20th percentile within the winning radius, because this round suggests the local niche is plausible while the exact tail definition may determine whether EC adds signal beyond the CA1 panel.
