## Summary
Tested the ca1-reactive-dominant-cuff family; candidate_variant_b won with selection score 0.5291.

## Metrics
Winner candidate_variant_b: partial r -0.5291, selection score 0.5291, LOO predictive r 0.4023, IS-LOO gap 0.1268, penalty 0.1268, adjusted score 0.4023.
Other tested variations: candidate_variant_a (selection 0.5204, partial r -0.5204, LOO 0.3780)

## Findings
1. What worked and why
   The winning feature measured the fraction of CA1 pyramidal neurons whose 35 um local astrocyte cuff had at least 2 astrocytes and met the reactive-dominance rule. The stricter supermajority threshold better isolates pyramidal neurons whose full local astrocyte cuff has shifted toward a reactive phenotype.
2. What failed and why
   The majority-reactive threshold looked too permissive: it counts mixed cuffs, so it likely dilutes the neuron-centered phenotype that stronger reactive dominance is isolating.
3. Error pattern: which donors are consistently wrong and what they share
   Largest absolute LOO errors were H20.33.018, H21.33.018, H21.33.023; 1/3 have Braak >=5. Example high-error donors: H20.33.018 (obs -0.30, pred -0.18), H21.33.018 (obs -0.14, pred -0.02), H21.33.023 (obs -0.15, pred -0.04)

## Rationale
This approach is biologically coherent because it stays within CA1, keeps the neuron as the center of the niche, and changes only the composition rule of the surrounding cuff.
Compared with nearby alternatives, the winner best balanced specificity against sparsity: it sharpened the cuff phenotype without making the donor score too unstable.
Because the current panel already contains broad astrocytic reactivity and count-based cuff features, this reactive-dominance score is most plausibly additive if cuff composition carries information beyond raw reactive-neighbor count; that final judgment belongs to the panel evaluator.

## Interpretation
The signal appears to reflect CA1 peripyramidal astrocyte niches in which the local astrocyte population has shifted toward reactive phenotype dominance around pyramidal neurons.
Population: CA1 pyramidal neurons with neighboring CA1 astrocytes/reactive astrocytes. Niche: 35 um peripyramidal cuff. Feature summary: donor-level fraction of neurons meeting the reactive-dominant cuff rule. Simplest observable pattern: CA1 pyramidal neurons wrapped by astrocyte cuffs in which reactive astrocytes clearly outnumber quiescent astrocytes.

## Next
Run the next local sweep by keeping the same CA1 peripyramidal niche but varying the minimum cuff size above 2 neighbors, since the current comparison mainly tests composition threshold and the remaining error may come from weakly populated cuffs.
