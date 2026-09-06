## Summary
One sentence: tested the ca1_reactive_neuron_exclusion family, reactive_pyramidal_desert_t20 won, and its local winner selection score was 0.2545.

## Metrics
Winning variation `reactive_pyramidal_desert_t20` used feature column `ca1_reactive_pyramidal_desert_t20` with IS partial r -0.2545, selection score 0.2545, LOO predictive r 0.0064, IS-LOO gap 0.2481, gap penalty 0.0490, and adjusted score -0.0426.

Other tested variations:
- reactive_pyramidal_desert_t20: partial_r=-0.2545 selection=0.2545 loo_r=0.0064
- reactive_pyramidal_desert_t10: partial_r=-0.1977 selection=0.1977 loo_r=-0.0625

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The winner, reactive_pyramidal_desert_t20, produced the strongest in-sample residual association in this local sweep by tracking the donor-level fraction of CA1 reactive astrocytes sitting in pyramidal-neuron-poor microenvironments within a 35 µm CA1 niche. That biological framing is coherent with a reactive glial response to local neuron depletion or scar-like remodeling.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The nearby alternative reactive_pyramidal_desert_t10 was weaker, and more importantly the family as a whole generalized poorly (winner LOO predictive r 0.0064). That suggests the thresholded desert fraction is biologically sensible but still too unstable or too redundant with existing CA1 injury features to yield reliable donor-level prediction on its own.
3. Error pattern: which donors are consistently wrong and what they share
- H20.33.018: outcome -0.304, prediction -0.106, ca1_reactive_pyramidal_desert_t20=0.663, Braak 6, CERAD 3, status NA.
- H21.33.040: outcome -0.033, prediction -0.215, ca1_reactive_pyramidal_desert_t20=0.670, Braak 5, CERAD 3, status NA.
- H21.33.023: outcome -0.155, prediction -0.044, ca1_reactive_pyramidal_desert_t20=0.649, Braak 4, CERAD 0, status NA.

## Rationale
reactive_pyramidal_desert_t20 won because the threshold of 20% is less extreme than the 10% cutoff, so it captures a broader band of neuron-poor reactive niches rather than only the most severe deserts. That likely explains why it outperformed the nearby alternative locally, although the weak LOO signal means the family still looks unstable or redundant.
Relative to the current panel, this niche-focused reactive-astrocyte measure is biologically more specific than a global reactive burden feature, but the near-zero LOO correlation makes it unlikely to add robust new information without another refinement.

## Interpretation
The signal most plausibly represents a CA1 reactive astrocyte niche in which local pyramidal neurons have dropped out, leaving reactive astrocytes embedded in neuron-poor tissue. Population: CA1 reactive astrocytes. Niche: their 35 µm CA1 neighborhood. Summary: fraction of reactive centers below a local pyramidal-neuron fraction threshold. Observable pattern: clusters of reactive astrocytes occupying locally neuron-sparse CA1 pockets.

## Next
Next local sweep: keep the same CA1 reactive-astrocyte-centered niche, but test whether weighting each center by its local reactive-cell density or requiring a minimum neighbor count sharpens the signal among the donors with the largest LOO residuals.
