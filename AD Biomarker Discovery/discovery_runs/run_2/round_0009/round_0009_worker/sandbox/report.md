## Summary
One sentence: tested the ca1_reactive_hotspot_corpora_amylacea_enrichment family, candidate_variant_a won, and its local winner selection score was 0.2571.

## Metrics
Winner: candidate_variant_a with partial r -0.2571, selection score 0.2571, LOO predictive r 0.1456, IS-LOO gap 0.4027, penalty 0.1007, and adjusted score 0.1564.
- candidate_variant_b: partial r -0.0363, selection score 0.0363, LOO r -0.1791.

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The working signal was the donor-level difference between CA1 reactive-astrocyte hotspot neighborhoods and matched astrocyte hotspot neighborhoods in local Corpora Amylacea fraction. This is biologically coherent because Corpora Amylacea can mark chronic astrocytic waste-handling or clearance-failure niches, so enrichment specifically around reactive hotspots should track more pathologic gliosis than bulk CA1 composition alone.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The nearby alternative candidate_variant_b at r=90 px was weaker (selection score 0.0363 vs 0.2571), suggesting that broadening the niche diluted the CA1 hotspot-local chronic-waste contrast.
3. Error pattern: which donors are consistently wrong and what they share
   - Largest absolute LOO errors: H20.33.018 (abs error 0.203, Braak 6, CERAD 3); H21.33.040 (abs error 0.119, Braak 5, CERAD 3); H21.33.023 (abs error 0.117, Braak 4, CERAD 0).
   - Shared pattern: Most high-error donors share Braak stage=5; Most high-error donors share CERAD=3.

## Rationale
The best variation kept the already productive CA1 hotspot framework from prior rounds but changed the neighboring population to Corpora Amylacea, making the scalar about chronic waste-accumulation content around reactive astrocyte microclusters rather than neuronal or oligodendroglial depletion. It beat the nearby alternative because its radius better matched the spatial scale at which reactive hotspots appear most different from homeostatic astrocyte hotspots in local chronic-debris burden. This looks plausibly additive to the current panel because it reads out a niche-content state rather than a direct abundance or morphology summary.

## Interpretation
The signal appears to mean that CA1 reactive astrocyte hotspots are embedded in Corpora-Amylacea-rich microenvironments relative to matched astrocyte hotspots. Population: Reactive Astrocyte versus Astrocyte hotspot centers. Niche: CA1 radius-70px local neighborhoods. Feature summary: mean smoothed local Corpora Amylacea fraction around reactive hotspots minus the same quantity around astrocyte hotspots. Simplest observable pattern: focal reactive astrocyte clusters sitting inside CA1 pockets with visibly more Corpora Amylacea.

## Next
Try a tighter CA1 corpora-amylacea niche sweep around the winning 70 px scale (for example 50-70 px, keeping k=3) to test whether the signal sharpens at a more focal hotspot radius.
