set -euo pipefail
cat > /scratch/report.md <<'MD'
## Summary
Tested the CA1 reactive-astrocyte plus lymphocyte niche family around CA1 pyramidal neurons; the local winner was `candidate_variant_b` (`ca1_reactive_lymphocyte_neighbor_fraction_r128px`) with selection score `0.3492`.

## Metrics
Winning variation: `candidate_variant_b` / `ca1_reactive_lymphocyte_neighbor_fraction_r128px`

- IS partial r: `-0.3492`
- Selection score: `0.3492`
- LOO predictive r: `-0.0169`
- IS-LOO Gap: `0.3323`
- Penalty: `0.3323`
- Adjusted Score: `0.0169`

Ranking of tested variations:
1. `candidate_variant_b` (`r128px`): partial r `-0.3492`, selection score `0.3492`, LOO predictive r `-0.0169`
2. `candidate_variant_a` (`r96px`): partial r `-0.3339`, selection score `0.3339`, LOO predictive r `-0.1106`

So the broader 128 px neighborhood slightly improved the in-sample niche signal over 96 px, but both variants had weak out-of-sample behavior and large IS-LOO gaps.

## Findings
1. What worked and why  
   The family did recover a moderate negative in-sample partial correlation: donors with a larger fraction of CA1 pyramidal neurons sitting in a local niche that contains both reactive astrocytes and lymphocytes tended to have faster memory decline. That is biologically coherent with an inflammatory subtype of the already-productive CA1 injury niche. The 128 px radius beat 96 px, suggesting the immune-associated context is slightly broader and sparser than the immediate reactive-pyramidal contacts already captured by the accepted panel.

2. What failed and why  
   The feature was too sparse and unstable to generalize. Sixteen of 35 donors had a zero-valued winning feature, and most nonzero values were still very small. Lymphocytes are rare in CA1, so this conjunction feature behaves like a rare-event filter layered on top of an already selective niche. That created moderate in-sample separation but almost no fold-stable predictive signal (`LOO r = -0.0169`). The broader radius helped coverage of the niche slightly, but not enough to overcome the sparsity problem.

3. Error pattern: which donors are consistently wrong and what they share  
   - **Undercalled severe-decline donors with little/no lymphocyte niche signal:** `H20.33.018`, `H20.33.038`, `H21.33.023`, `H20.33.037`, `H20.33.046` all had substantially worse memory-decline outcomes than predicted. Several have zero or near-zero winning feature values despite marked decline, implying the damaging CA1 state in these donors is not dependent on explicit local lymphocyte presence.
   - **Overcalled donors with lymphocyte-rich local signal but milder decline:** `H20.33.034`, `H20.33.013`, `H21.33.040`, `H21.33.035` were predicted worse than observed. In particular `H20.33.034` had unusually high CA1 lymphocyte count and a high niche fraction, but only modest decline, suggesting that diffuse or incidental lymphocyte presence can inflate the biomarker without tracking the core neuron-loss program.

## Rationale
The best variation is biologically coherent because it asks whether the established CA1 reactive astrocyte–pyramidal injury niche becomes more pathologic when an immune-cell component is also nearby. It beat the 96 px alternative because the relevant lymphocyte context appears to operate on a slightly broader spatial scale than the tighter reactive-neighbor definitions used in earlier rounds.  

However, this family does **not** currently look likely to add useful new information beyond the accepted panel. Its in-sample direction is sensible, but the almost-null LOO signal suggests that explicit lymphocyte conditioning is mostly slicing the existing CA1 injury niche into a sparse subset rather than identifying a robust new donor-level axis.

## Interpretation
The signal seems to mean: **a subset of CA1 pyramidal neurons sit in inflammatory reactive-astrocyte niches that also contain nearby lymphocytes, and those donors tend to have faster decline in-sample**.  

- **Population:** CA1 pyramidal neurons, CA1 reactive astrocytes, CA1 lymphocytes  
- **Niche:** CA1 pyramidal-centered neighborhood requiring both reactive astrocyte and lymphocyte presence  
- **Feature summary:** donor-level fraction of CA1 pyramidal neurons meeting that dual-neighbor condition  
- **Simplest observable pattern:** scattered CA1 pyramidal zones where reactive gliosis coincides with focal lymphocyte presence  

Biologically, that looks like an inflammatory subtype of the CA1 reactive injury niche, but in this cohort it is too rare and inconsistent to serve as a stable standalone biomarker.

## Next
Because `r128px` beat `r96px` only slightly and both variants were dominated by sparsity, the next local sweep should keep the broader CA1 immune context but change the donor summary away from a rare pyramidal fraction—for example, test a **CA1 reactive-astrocyte–anchored lymphocyte co-localization rate at 128 px** to see whether inflammatory burden is better captured on the glial side rather than by counting the small subset of affected pyramidal anchors.
MD
echo 'wrote /scratch/report.md'
