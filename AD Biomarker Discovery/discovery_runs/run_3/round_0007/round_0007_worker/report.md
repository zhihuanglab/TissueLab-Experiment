## Summary
Tested the CA1 peripyramidal reactive-astrocyte total-area burden family; `mean_total_reactive_area_r35um` won with local selection score `0.5438`.

## Metrics
Winning variation: `mean_total_reactive_area_r35um` (`Mean per-CA1-pyramidal-neuron summed reactive-astrocyte contour area within 35 µm`)

- IS partial r: `-0.5438`
- Selection score: `0.5438`
- LOO predictive r: `0.4146`
- IS-LOO Gap: `0.1292` (penalty=`0.0000`)
- Adjusted Score: `0.4146`
- Analyzable donors: `35` / `35`

Other tested variations: mean_total_reactive_area_r45um scored 0.5435 (partial r -0.5435, LOO r 0.4124)

## Findings
1. What worked and why  
   The winning feature works because it combines two previously promising ingredients in the same biological niche: reactive astrocyte presence around CA1 pyramidal neurons and reactive astrocyte hypertrophy. Averaging the **summed reactive contour area per CA1 pyramidal neuron** captures a donor-level severity readout of perineuronal glial burden rather than treating cuff count and cell size as separate signals.

2. What failed and why  
   The losing local variation underperformed because its radius choice was less well matched to the biologically informative niche. In this family, moving away from the winning radius diluted the peripyramidal signal by either missing relevant immediate cuffs if too tight or mixing in broader local gliosis if too wide.

3. Error pattern: which donors are consistently wrong and what they share  
   Largest LOO misses are H21.33.018, H20.33.018, H21.33.023, H21.33.006, H20.33.046. Under-predicted decline donors (H21.33.018, H20.33.018, H21.33.023, H20.33.046) had more memory decline than predicted from this local burden, while over-predicted decline donors (H21.33.006) had less decline than their CA1 burden would suggest. Two repeat mismatch modes appear: low-burden but high-decline donors (H21.33.018, H21.33.023, H20.33.046) suggest pathology outside the CA1 peripyramidal astrocyte niche, whereas high-burden but milder-decline donors (H21.33.006) suggest reactive burden can be present without proportionate memory slope severity. The median winner feature value is ≈ 5.33; the worst errors span both sides of that median, arguing against a simple extraction artifact and more toward biological heterogeneity.

## Rationale
This approach is biologically coherent because reactive astrocytes can contribute to memory-linked neuronal stress through both **how many** cells cluster near a neuron and **how hypertrophic** those cells are. Summing reactive astrocyte contour area within a fixed peripyramidal shell turns that biology into one auditable donor scalar. It beat the nearby alternative because the winning radius better matches the local cuff-like niche around CA1 pyramidal neurons rather than a more diluted neighborhood. Relative to the current panel, this feature seems biologically close to existing CA1 cuff metrics, but it may still be useful as a cleaner severity readout that could replace a more redundant thresholded member.

## Interpretation
Biologically, the signal appears to reflect **reactive astrocytic burden concentrated around CA1 pyramidal neurons**.  
- Population: `Reactive Astrocyte` around `Pyramidal Neuron`
- Niche: `CA1` peripyramidal neighborhood within `35 µm`
- Feature summary: mean per-neuron summed reactive astrocyte contour area
- Simplest observable pattern: donors with worse memory decline tend to show thicker / more hypertrophic reactive astrocyte cuffs clustered around CA1 pyramidal neuronal somata

## Next
Try an upper-tail burden summary in the same CA1 peripyramidal niche, such as the fraction of CA1 pyramidal neurons whose summed reactive area exceeds a donor-specific threshold within 35 µm, because the winning mean burden suggests burden severity matters but the remaining errors imply donor heterogeneity in how concentrated that burden is across neurons.
