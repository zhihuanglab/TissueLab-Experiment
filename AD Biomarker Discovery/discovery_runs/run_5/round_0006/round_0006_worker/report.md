## Summary
One sentence: tested the CA1 reactive-astro-associated pyramidal isolation family; candidate_variant_a won, with selection score 0.5244.

## Metrics
Winning variation `candidate_variant_a` (`ca1_pyramidal_reactive_astro_isolated_fraction_k5_30um`) had partial r -0.5244, selection score 0.5244, LOO predictive r 0.3922, IS-LOO gap 0.1322, penalty 0.1322, and adjusted score 0.3922. Other tested variations ranked as: candidate_variant_a: partial_r=-0.5244, selection_score=0.5244, loo_predictive_r=0.3922; candidate_variant_b: partial_r=-0.5200, selection_score=0.5200, loo_predictive_r=0.3804.

## Findings
1. What worked and why (tie to the biological meaning of the target): The signal came from **CA1 pyramidal neurons** and sharpened when restricted to those with a **reactive astrocyte within 30 um** while also having **few nearby CA1 pyramidal peers**. That donor-level fraction appears to summarize local neuronal depletion inside the same reactive-astro niche highlighted by earlier accepted features.
2. What failed and why (specific to the chosen hypothesis and what went wrong): candidate_variant_b underperformed because relaxing the isolation threshold to eight nearby pyramidal neighbors diluted the signal with less-depleted CA1 neighborhoods.
3. Error pattern: Largest LOO misses were H20.33.018 (abs err 0.141, braak 6, cerad 3); H21.33.018 (abs err 0.139, braak 2, cerad 3); H21.33.023 (abs err 0.128, braak 4, cerad 0); their pathology burden is mixed rather than uniform (braak range 2-6, cerad range 0-3), suggesting this niche is not just a proxy for global AD stage.

## Rationale
The k<=5 rule keeps the biomarker focused on genuinely sparse CA1 pyramidal neighborhoods, which is the closest match to reactive-astro-associated local dropout. It beat the k<=8 alternative because the looser cutoff admits neurons that still have several nearby pyramidal peers, so it is less specific to tissue disorganization.
If this candidate is kept, it is likely to add more neuron-state specificity than a pure astrocyte abundance summary because it encodes **which CA1 pyramidal neurons sit inside the reactive-astro niche and how locally isolated they are**.

## Interpretation
Biologically, the signal seems to mean **reactive-astro-associated local CA1 pyramidal dropout / sparsification**.  
Population: **CA1 Pyramidal Neuron**.  
Niche: **within 30 um of a Reactive Astrocyte in CA1**.  
Feature summary: **fraction of CA1 pyramidal neurons meeting the reactive-astro-near and low-neighbor criterion**.  
Simplest observable pattern: **more lone or sparsely surrounded CA1 pyramidal somata sitting next to reactive astrocytes**.

## Next
Tighten the same family around the winning side of the sweep: test k=3-6 or require at least two reactive astrocytes within 30 um, because the looser k=8 cutoff lost and the remaining errors likely reflect donors where astrocyte exposure is present but neuronal depletion must be defined more sharply.
