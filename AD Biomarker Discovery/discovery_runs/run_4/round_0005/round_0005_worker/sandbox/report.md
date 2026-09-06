## Summary
Tested the CA1 reactive-astrocyte local purity family; the 60 µm neighborhood variation (`reactive_purity_r60um`) won with selection score 0.4781.

## Metrics
Winning variation: `reactive_purity_r60um` (`feature_column = ca1_reactive_local_purity_r60um`).

- IS partial r: -0.4781
- Selection score: 0.4781
- LOO predictive r: 0.3323
- IS-LOO Gap: 0.1458
- Penalty: 0.0000
- Adjusted Score: 0.3323
- p-value: 0.0037
- n analyzable / total: 35 / 35

Ranking of tested variations:
1. `reactive_purity_r60um` — partial r -0.4781, selection 0.4781, LOO r 0.3323
2. `reactive_purity_r40um` — partial r -0.4613, selection 0.4613, LOO r 0.3199

The broader 60 µm patch scale beat the 40 µm baseline on both in-sample partial correlation and LOO predictive correlation, while keeping the IS-LOO gap below the 0.15 penalty threshold.

## Findings
1. What worked and why  
   The strongest signal came from **CA1 Reactive Astrocytes specifically**, conditioned on their **local astrocyte-lineage niche** rather than on global abundance alone. Donors with higher mean local reactive purity around CA1 Reactive Astrocytes had faster memory decline after adjusting for age, Braak, CERAD, and sex. The 60 µm radius likely captures a biologically coherent **reactive glial patch** scale better than the tighter 40 µm radius.

2. What failed and why  
   The 40 µm version was slightly weaker. It appears too local and noisier: it produced fewer eligible centers on average (`~185` valid reactive centers per donor vs `~417` at 60 µm), so the donor summary is less stable. That suggests the relevant lesion architecture is not just immediate cell-to-cell adjacency, but a somewhat broader reactive microdomain.

3. Error pattern: which donors are consistently wrong and what they share  
   Largest LOO misses included `H20.33.018`, `H21.33.023`, `H21.33.006`, and `H21.33.018`. Two opposite error modes appear:
   - **Very high local purity but decline not as extreme as predicted**: `H21.33.006` (purity 0.8402, predicted much worse than observed).
   - **Lower local purity but decline worse than predicted**: `H21.33.023` (purity 0.3850) and `H21.33.018` (purity 0.3198).
   
   So this biomarker captures one injury architecture well, but some donors likely decline through a different mechanism not dominated by reactive-astrocyte patch purity. The worst underpredicted severe donors (`H20.33.018`, `H20.33.046`, `H20.33.037`) still have moderate-to-high purity, implying local purity is relevant but not sufficient to explain the full severity range.

## Rationale
This winner is biologically coherent because it tests whether CA1 reactive astrocytes occur as **locally self-reinforcing reactive islands** instead of being diffusely intermixed with non-reactive astrocytes. That is a plausible tissue correlate of focal glial injury architecture. It beat the nearby 40 µm alternative because the broader neighborhood better captures contiguous astrocyte-lineage reactivity and yields more stable donor summaries. Relative to the current panel, it looks plausibly additive because it is not just another global reactive-astrocyte count: it measures **how reactive those astrocytes are spatially organized within CA1 astrocyte neighborhoods**.

## Interpretation
The signal seems to mean that **memory decline tracks CA1 reactive astrocytes that cluster into locally reactive astrocyte-lineage patches**.  

- **Population:** CA1 Reactive Astrocytes, with Astrocytes + Reactive Astrocytes used as the local denominator  
- **Niche:** CA1 local astrocyte-lineage neighborhood within ~60 µm  
- **Feature summary:** donor-level mean local reactive purity around each CA1 Reactive Astrocyte, requiring at least 2 local astrocyte-lineage neighbors  
- **Observable tissue pattern:** patches where a reactive astrocyte is surrounded mostly by other reactive astrocytes rather than by mixed reactive/non-reactive astrocytes

## Next
Do one more local sweep inside this same family by keeping the 60 µm neighborhood and testing a **high-purity burden** summary, such as the fraction of CA1 Reactive Astrocytes whose local reactive purity exceeds 0.7, because 60 µm beat 40 µm and the main residual errors suggest that mean purity may miss donors with a smaller number of especially pure reactive patches.
