## Summary
Tested the CA1-region pyramidal-neuron fraction family with a CA1 vs CA1+CA2 localization sweep; `candidate_variant_a` (CA1 only) won with selection score 0.2940.

## Metrics
Winning variation: `candidate_variant_a` / `ca1_pyramidal_fraction`.

- IS partial r: 0.2940
- Selection score: 0.2940
- LOO predictive r: 0.0976
- IS-LOO gap: 0.1964
- Penalty: 0.1964
- Adjusted score: 0.0976
- Coverage: 35/35 donors
- p-value: 0.0865

Ranking of nearby variations:
1. `candidate_variant_a` (CA1 only): partial r 0.2940, selection score 0.2940, LOO predictive r 0.0976
2. `candidate_variant_b` (CA1+CA2 union): partial r 0.2668, selection score 0.2668, LOO predictive r 0.0573

So the local sweep favored tighter CA1 localization; adding CA2 reduced both in-sample and LOO signal.

## Findings
1. What worked and why  
   The simple donor-level fraction of classified CA1 cells labeled `Pyramidal Neuron` gave the stronger local signal. This is biologically aligned with the memory target because CA1 pyramidal neurons are a core hippocampal output population, so a lower CA1 neuronal share plausibly reflects local neuronal loss or local replacement by glial / non-neuronal classes in a memory-relevant niche.

2. What failed and why  
   Broadening the niche from CA1 to CA1+CA2 weakened the signal. That suggests the effect is more CA1-specific than a generic cornu ammonis neuronal-composition change. The family also shows limited predictive robustness: the winner's LOO r (0.0976) is much smaller than its in-sample partial r (0.2940), so this scalar captures a real directional trend but not the full donor-to-donor heterogeneity in memory decline.

3. Error pattern: which donors are consistently wrong and what they share  
   The same donors are among the largest errors for both CA1-only and CA1+CA2:
   - `H20.33.018` is the main underpredicted severe-decline donor in both variants and is also the only unstable donor by leave-one-out influence.
   - `H21.33.023`, `H21.33.004`, and `H21.33.046` are overpredicted as relatively preserved because they retain high CA1 pyramidal fractions despite worse-than-expected memory decline.
   - `H21.33.040` goes the other direction: middling CA1 fraction but much better outcome than predicted.  
   A visible pattern is that several large-error donors with unexpectedly poor outcome despite high CA1 pyramidal fraction are male donors, implying this biomarker misses decline mechanisms that do not primarily express as CA1 neuronal depletion.

## Rationale
The best variation is biologically coherent because it measures the representation of the canonical excitatory neuron class inside the hippocampal subfield most closely tied to episodic memory vulnerability. It beat the CA1+CA2 alternative because CA2 likely dilutes the signal with a neighboring region whose cellular vulnerability profile is less tightly coupled to the memory phenotype. With the panel currently empty, this is a plausible seed candidate because it is simple, auditable, high-coverage, and directly interpretable, although its modest LOO performance means it should be viewed as a tentative anchor rather than a stable standalone predictor.

## Interpretation
The signal appears to mean: donors with a more neuron-sparse CA1 compartment tend to have worse memory slope after adjustment for age, sex, Braak stage, and CERAD burden.  

- Population: `Pyramidal Neuron`
- Niche: CA1 annotated hippocampal region
- Feature summary: donor-level fraction `n(CA1 pyramidal neurons) / n(all classified CA1 cells)`
- Simplest observable pattern: CA1 looks relatively depleted of pyramidal somata and correspondingly more occupied by non-pyramidal cells on slides from donors with worse memory decline.

Because `slope_zmem0` is negative-valued in this cohort, the positive partial correlation implies the converse of the original question is what the data support locally: higher CA1 pyramidal fraction tracks relatively better memory slope, and lower fraction tracks worse decline.

## Next
Stay in the CA1 pyramidal-fraction family, but test denominator cleanup rather than region expansion: compare the CA1 pyramidal fraction using all classified cells versus excluding `Negative control` class from the denominator, because CA1 beat CA1+CA2 and the remaining error pattern suggests denominator contamination / non-biological class mixture may be weakening robustness.
