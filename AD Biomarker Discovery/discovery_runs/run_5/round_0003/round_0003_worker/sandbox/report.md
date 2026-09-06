## Summary
Tested the **CA1 peripyramidal reactive-astrocyte crowding** family; the winner was **k3_within_30um** with local winner selection score **0.5247**.

## Metrics
Winning variation: **k3_within_30um**  
- IS partial r: **-0.5247**
- Selection score: **0.5247**
- LOO predictive r: **0.4026**
- IS-LOO Gap: **0.1220** with penalty **0.1220**
- Adjusted Score: **0.4026**

Local ranking:
1. **k3_within_30um** — partial r **-0.5247**, selection score **0.5247**, LOO predictive r **0.4026**
2. **k2_within_30um** — partial r **-0.5084**, selection score **0.5084**, LOO predictive r **0.3492**

So the stricter threshold (at least 3 nearby CA1 reactive astrocytes within ~30 um of a CA1 pyramidal neuron) beat the milder crowding threshold (at least 2 neighbors) on both in-sample and LOO diagnostics.

## Findings
1. **What worked and why**  
   The signal strengthened when the niche was sharpened to **rare, high-density reactive astrocyte crowding directly around CA1 pyramidal neurons**. Requiring **3+** nearby reactive astrocytes likely enriches for a more pathologic perineuronal glial state than the broader k2 definition. The negative partial correlation means donors with a larger fraction of these crowded pyramidal-neuron niches tend to have worse memory-decline outcome values.

2. **What failed and why**  
   The weaker **k2** threshold still worked, but it diluted the signal by counting milder exposure states that are probably closer to background reactive-astro presence. This family also remains inherently close to the already accepted round-2 marker (any reactive-astro neighbor within 30 um), so the main risk is **redundancy rather than lack of signal**.

3. **Error pattern: which donors are consistently wrong and what they share**  
   Largest LOO errors came from donors such as **H20.33.018, H20.33.037, H21.33.043, H21.33.023, and H20.33.046**. Two patterns stood out:
   - **Underpredicted severe decline with little k3 crowding**: H20.33.037, H21.33.023, and H20.33.046 had fairly negative outcomes despite near-zero k3 feature values, suggesting additional non-astrocytic or non-peripyramidal injury mechanisms.
   - **Overcalled decline from unusually high crowding**: H21.33.043 had one of the highest k3 crowding values but only modest decline, suggesting that this niche is not sufficient on its own and may need a morphology or cell-state refinement.
   - **Magnitude miss at the extreme high end**: H20.33.018 had the highest k3 value and strong decline in the expected direction, but the model still underpredicted how severe the decline was.

## Rationale
This winner is biologically coherent because it focuses on a very specific microenvironment: **CA1 pyramidal neurons embedded in dense clusters of nearby reactive astrocytes**. That is a sharper version of the accepted reactive-astro lead than either overall CA1 reactive-astro burden or binary any-neighbor exposure. It beat k2 because the stricter threshold appears to isolate a more severe local glial response rather than merely detecting incidental proximity.

Because the accepted panel already contains both **CA1 reactive astrocyte fraction** and **CA1 pyramidal any-neighbor reactive-astro exposure**, this round's winner seems most likely to contribute as a **refinement or possible replacement candidate**, not obviously as a strongly orthogonal third signal. Still, its slightly stronger single-feature partial r and better LOO diagnostic suggest that severe crowding is the cleaner formulation within this local family.

## Interpretation
The signal appears to reflect **focal reactive-astrocyte encirclement of CA1 pyramidal neurons**.

- **Population:** CA1 Pyramidal Neurons with CA1 Reactive Astrocytes as neighbors
- **Niche:** Peripyramidal CA1 microenvironment within ~30 um
- **Feature summary:** Fraction of CA1 pyramidal neurons having **at least 3** nearby CA1 reactive astrocyte centroids
- **Simplest observable tissue pattern:** rare CA1 neuronal neighborhoods where a pyramidal soma is surrounded by a small cluster of reactive astrocytes rather than a single adjacent glial cell

Biologically, that pattern is consistent with a more intense local glial reactivity state that tracks worse memory decline.

## Next
Run one nearby sweep that replaces the hard threshold with a **continuous per-pyramidal reactive-astro burden score** in the same CA1 niche, for example the **mean number of CA1 reactive astrocytes within 30 um per CA1 pyramidal neuron** or an upper-tail summary. Since **k3 beat k2**, severity seems to matter; a continuous burden measure may preserve that severity information while better handling donors like H20.33.018 (extreme burden) and separating them from modest-high but less consequential cases like H21.33.043.
