## Summary
One sentence: tested the ca1 pyramidal neuron abundance family, `ca1_pyramidal_fraction` won, and its local winner selection score was 0.2940.

## Metrics
Winning variation `ca1_pyramidal_fraction` had partial r 0.2940, selection score 0.2940, LOO predictive r 0.1999, IS-LOO gap 0.0940, penalty 0.0940, and adjusted score 0.1999.
Other tested variations:
- `ca1_pyramidal_density`: partial r 0.2219, selection score 0.2219, LOO predictive r 0.1856.

## Findings
1. What worked and why (tie to the biological meaning of the target)
   - The winner tracked donor-to-donor variation in CA1 pyramidal abundance while controlling for age, Braak stage, CERAD, and sex. This is biologically aligned with memory-circuit vulnerability because CA1 pyramidal neurons are a direct substrate of hippocampal memory output.
2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The losing normalization was more sensitive to regional annotation size, so it mixed neuronal abundance with donor-to-donor CA1 area variability. That diluted the cell-loss signal relative to the cell-composition fraction.
3. Error pattern: which donors are consistently wrong and what they share
   - Largest absolute LOO errors were: H20.33.018 (outcome -0.3042, predicted -0.1133, ca1_pyramidal_fraction 0.1379); H21.33.040 (outcome -0.0328, predicted -0.1696, ca1_pyramidal_fraction 0.1931); H21.33.023 (outcome -0.1546, predicted -0.0209, ca1_pyramidal_fraction 0.2494). These donors likely share mismatch between CA1 pyramidal abundance and decline severity, suggesting additional gliosis, extra-CA1 pathology, or circuit-level compensation not captured by a single abundance scalar.

## Rationale
The best variation is biologically coherent because it isolates the vulnerable CA1 pyramidal population inside its native hippocampal niche and expresses it as a simple donor-level scalar. It beat the nearby density alternative because dividing by total CA1 cells better measures selective compositional loss, whereas density is partially confounded by how much CA1 was annotated and how expanded or contracted the region appears on the section. As a seed biomarker, it is likely to add interpretable information because the current accepted panel is empty.

## Interpretation
The signal seems to mean biologically that donors with lower CA1 pyramidal abundance tend to have faster memory decline. The population is CA1 Pyramidal Neuron cells, the niche is annotated CA1 pyramidal layer territory within the hippocampus, the feature summary is donor-level fraction of CA1 cells classified as Pyramidal Neuron, and the simplest observable pattern is a CA1 field with visibly fewer pyramidal neurons relative to the other annotated CA1 cellular content. Lowest-feature donors: H21.33.012 (ca1_pyramidal_fraction 0.0862), H20.33.038 (ca1_pyramidal_fraction 0.1048), H21.33.043 (ca1_pyramidal_fraction 0.1160). Highest-feature donors: H21.33.046 (ca1_pyramidal_fraction 0.3148), H21.33.033 (ca1_pyramidal_fraction 0.2826), H20.33.024 (ca1_pyramidal_fraction 0.2790).

## Next
One specific suggestion for the next local sweep based on the error pattern and which nearby variations won or lost: stay in CA1 and test whether reactive astrocyte enrichment or astrocyte-to-pyramidal balance within the same CA1 niche explains the donors that remain badly predicted after the pyramidal-fraction winner.
