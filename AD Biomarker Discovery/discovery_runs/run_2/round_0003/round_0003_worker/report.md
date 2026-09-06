## Summary
One sentence: tested the ca1_pyramidal_reactive_niche family, the winning local variation was ca1_pyr_reactive_niche_r70px, and its local winner selection score was 0.4837.

## Metrics
Winning variation: ca1_pyr_reactive_niche_r70px
- IS partial r: -0.4837
- Selection score: 0.4837
- LOO predictive r: 0.3360
- IS-LOO gap: 0.1477
- Gap penalty: 0.0000
- Adjusted score: 0.3360
Other tested radii: ca1_pyr_reactive_niche_r60px (selection=0.4745, partial_r=-0.4745, loo=0.3217); ca1_pyr_reactive_niche_r50px (selection=0.4643, partial_r=-0.4643, loo=0.3059).

## Findings
1. What worked and why: ca1_pyr_reactive_niche_r70px was the strongest radius, suggesting the signal is most coherent when the donor scalar measures the fraction of CA1 pyramidal neurons whose astroglial-contacting neighborhood specifically contains a reactive astrocyte. In analyzable donors the mean feature value was 0.5333, with mean reactive-neighbored pyramidal count 360.1 over mean astroglial-neighbored pyramidal count 723.3; this keeps the score anchored to a neuron-adjacent gliosis niche rather than bulk gliosis alone.
2. What failed and why: the nearby radius alternatives underperformed (ca1_pyr_reactive_niche_r60px (too tight), ca1_pyr_reactive_niche_r50px (too tight)). The tighter radius likely misses reactive astrocytes sitting just outside the immediate soma-scale shell, while the broader radius dilutes the neuron-adjacent niche into more generic local astroglial abundance.
3. Error pattern: The largest errors did not map cleanly to a single available metadata label. Top absolute errors: H20.33.018 (abs err=0.165, status=unknown, Braak=6, CERAD=3); H21.33.023 (abs err=0.158, status=unknown, Braak=4, CERAD=0); H21.33.018 (abs err=0.132, status=unknown, Braak=2, CERAD=3).

## Rationale
This winner is biologically coherent because it conditions reactive astrocyte presence on surviving CA1 pyramidal neurons and normalizes by local astroglial availability, so it asks whether gliosis is preferentially concentrated right around the neuronal population already implicated by earlier rounds. It beat the nearby alternatives because the winning radius appears to balance specificity and capture: small enough to stay neuron-adjacent, but large enough to include the local reactive astrocyte shell. It is likely at least partially distinct from a pure bulk reactive fraction because the denominator removes non-specific astroglial abundance and centers the summary on pyramidal neighborhoods.

## Interpretation
The signal seems to mean a CA1 neuron-adjacent reactive gliosis state: among CA1 pyramidal neurons that sit in an astroglial neighborhood, donors with lower values of the reactive-over-astroglial neighborhood fraction showed lower niche scores tracked lower slope_zmem0. Population: CA1 Pyramidal Neuron query cells with CA1 Reactive Astrocyte / Astrocyte neighbors. Niche: the immediate peri-neuronal CA1 astroglial neighborhood. Feature summary: reactive-neighbored pyramidal neurons divided by astroglial-neighbored pyramidal neurons. Simplest observable pattern: more of the surviving CA1 pyramidal neurons are ringed by reactive astrocytes rather than only generic astrocytes.

## Next
Next local sweep: keep the same CA1 pyramidal-centered niche but test a slightly finer radius bracket around the winner (for example ±10 px around 70px) or replace the binary neighbor hit with a nearest-reactive distance / local reactive count summary, because the current ranking suggests spatial scale matters more than the broader family definition.
