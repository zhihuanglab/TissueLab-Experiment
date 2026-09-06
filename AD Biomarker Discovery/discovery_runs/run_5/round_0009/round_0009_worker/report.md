## Summary
Tested the CA1 reactive-astro-exposed pyramidal lymphocyte-enrichment family; `candidate_variant_a` won, with local selection score 0.1593.

## Metrics
Winning variation: `candidate_variant_a` (`ca1_ra_exposed_pyramidal_lym_ge1_fraction_25um`).

- IS partial r: 0.1593
- Selection score: 0.1593
- LOO predictive r: 0.0657
- IS-LOO Gap: 0.0936
- Penalty: 0.0936
- Adjusted Score: 0.0657

Ranking of tested variations:
1. `candidate_variant_a` — partial r 0.1593, selection score 0.1593, LOO predictive r 0.0657
2. `candidate_variant_b` — all metrics `NaN` because the feature collapsed to a zero-variance donor vector

Coverage was effectively full for the winning feature, but the signal was very sparse: only 6 of 35 donors had a nonzero value.

## Findings
1. **What worked and why**
   - Conditioning on the **CA1 pyramidal neurons already within 25 um of a reactive astrocyte** did produce a weak positive signal when the readout was the **fraction of that exposed subset with at least one nearby lymphocyte**.
   - Biologically, this keeps the population and niche aligned with the previous panel direction: pyramidal neurons in a reactive-astrocyte microenvironment, sharpened further by local immune presence.
   - The donor scalar that carried the signal was `ca1_ra_exposed_pyramidal_lym_ge1_fraction_25um`, i.e. the fraction of reactive-astro-exposed CA1 pyramidal neurons that also had at least one lymphocyte within the same 25 um radius.

2. **What failed and why**
   - The stricter `>=2 lymphocytes` threshold failed completely. It produced a zero-variance feature across the cohort, so its partial correlation and LOO diagnostics were undefined.
   - This shows that the triple-conditioned niche is already rare at 25 um, and requiring two lymphocytes is too strict for this cohort.
   - Even the winning `>=1` feature remained extremely sparse: median 0, upper quartile 0, max only 0.0099, with nonzero values in just 6 donors.

3. **Error pattern**
   - The largest LOO misses were donors such as `H20.33.018`, `H21.33.023`, `H20.33.037`, and `H20.33.038`, all of which had **feature value 0** despite relatively poor memory-decline outcome.
   - The opposite error also appeared in milder-outcome donors with zero feature values, e.g. `H21.33.040`, `H21.33.035`, and `H20.33.024`, which were predicted more negative than observed.
   - The common pattern is that once the feature is zero, prediction falls back to confounds plus a very small residual correction. So most large errors are donors with **no detectable lymphocyte-positive triad at this strict niche definition**, meaning the biomarker under-covers the broader biology driving decline.

## Rationale
The best variation is biologically coherent because it isolates a specific inflammatory niche: **CA1 pyramidal neurons near reactive astrocytes and simultaneously near lymphocytes**. That is a plausible readout of focal neuroinflammatory stress in a vulnerable hippocampal compartment.

It beat the nearby alternative because `>=1 lymphocyte` preserved at least some donor-to-donor dynamic range, while `>=2 lymphocytes` was too rare to measure. The result suggests the niche is real but too sparse when quantified as a hard higher-count threshold. Given that rounds 7 and 8 already captured related immune-contact biology, this winner looks more like a **sharpened but weakly expressed refinement** than a clearly stronger new source of panel information.

## Interpretation
The signal appears to reflect **rare CA1 peripyramidal inflammatory foci** where a pyramidal neuron sits within 25 um of a reactive astrocyte and also has a nearby lymphocyte.

- **Population:** CA1 Pyramidal Neuron
- **Niche:** reactive-astrocyte-exposed, lymphocyte-positive local microenvironment within 25 um
- **Feature summary:** fraction of exposed pyramidal neurons with at least 1 nearby lymphocyte
- **Simplest observable tissue pattern:** small CA1 pockets where reactive astrocytes cluster around pyramidal neurons and occasional lymphocytes enter the same immediate neighborhood

Biologically, that likely marks a focal immune-reactive niche rather than a broad tissue state. The weakness of the donor-level signal suggests these foci may matter, but they are too infrequent in this exact binary definition to explain much additional variance on their own.

## Next
Keep the same CA1 reactive-astro-exposed pyramidal denominator, but next sweep a **continuous mean lymphocyte count within 25–30 um per exposed pyramidal neuron** rather than another hard threshold; `>=1` retained slight signal, while `>=2` failed by collapsing to zero.
