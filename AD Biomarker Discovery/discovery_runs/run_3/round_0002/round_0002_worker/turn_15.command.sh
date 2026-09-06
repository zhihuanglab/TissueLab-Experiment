python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
old = """    report = f\"\"\"## Summary
One sentence: tested the {HYPOTHESIS_FAMILY} family, {best['name']} won, and its local winner selection score was {best['selection_score']:.4f}.

## Metrics
Winning variation `{best['name']}` used feature column `{best['feature_column']}` with partial r {best['partial_r']:.4f}, selection score {best['selection_score']:.4f}, LOO predictive r {best['loo_predictive_r']:.4f}, IS-LOO gap {best['is_loo_gap']:.4f}, penalty {best['gap_penalty']:.4f}, and adjusted score {best['adjusted_score']:.4f}.
{\"\".join(line + chr(10) for line in other_lines)}

## Findings
1. What worked and why (tie to the biological meaning of the target)  
   The winner sharpened the accepted astrocytic-reactivity biology to a CA1 neuron-adjacent niche: reactive astrocytes among CA1 astrocyte-lineage cells lying within {best['radius_um']:.0f} µm of a CA1 pyramidal neuron. This directly targets peripyramidal gliosis around the neuron population most plausibly linked to memory decline.
2. What failed and why (specific to the chosen hypothesis and what went wrong)  
   The nearby radius alternative underperformed, suggesting the signal is not improved by simply broadening the neighborhood. If the broader radius lost, the niche is probably fairly tight; if the tighter radius lost, the signal likely needs a slightly more permissive peripyramidal field.
3. Error pattern: which donors are consistently wrong and what they share  
{summary_text}
{\"\".join(line + chr(10) for line in top_error_lines)}

## Rationale
The best variation is biologically coherent because it asks whether the CA1 astrocyte-reactivity signal is concentrated specifically around CA1 pyramidal neurons rather than diffusely across all of CA1. That beat the nearby alternative because the winning radius better matched the scale at which astrocyte activation appears to covary with donor memory-decline burden. Since the current accepted panel already contains a region-wide astrocytic reactivity member, this niche-refined candidate is likely partly redundant but could still add value if the evaluator finds that neuron-adjacent gliosis carries cleaner or slightly more specific signal.

## Interpretation
The signal seems to represent CA1 peripyramidal gliosis: the population is CA1 astrocyte-lineage cells, the niche is immediate proximity to CA1 pyramidal neurons, the donor-level summary is the reactive fraction among those niche astrocytes, and the simplest observable pattern is more reactive astrocytes crowding the CA1 pyramidal field in donors with worse memory decline.

## Next
Run one more local sweep around the winning niche scale by testing a slightly tighter or slightly broader CA1 pyramidal-neighborhood threshold and, if needed, a denominator-stabilized variant such as a shrinkage-adjusted reactive fraction.
\"\"\"
"""
new = """    metrics_lines = "\\n".join(other_lines)
    error_lines = "\\n".join(top_error_lines)
    report = f\"\"\"## Summary
One sentence: tested the {HYPOTHESIS_FAMILY} family, {best['name']} won, and its local winner selection score was {best['selection_score']:.4f}.

## Metrics
Winning variation `{best['name']}` used feature column `{best['feature_column']}` with partial r {best['partial_r']:.4f}, selection score {best['selection_score']:.4f}, LOO predictive r {best['loo_predictive_r']:.4f}, IS-LOO gap {best['is_loo_gap']:.4f}, penalty {best['gap_penalty']:.4f}, and adjusted score {best['adjusted_score']:.4f}.
{metrics_lines}

## Findings
1. What worked and why (tie to the biological meaning of the target)  
   The winner sharpened the accepted astrocytic-reactivity biology to a CA1 neuron-adjacent niche: reactive astrocytes among CA1 astrocyte-lineage cells lying within {best['radius_um']:.0f} µm of a CA1 pyramidal neuron. This directly targets peripyramidal gliosis around the neuron population most plausibly linked to memory decline.
2. What failed and why (specific to the chosen hypothesis and what went wrong)  
   The nearby radius alternative underperformed, suggesting the signal is not improved by simply broadening the neighborhood. If the broader radius lost, the niche is probably fairly tight; if the tighter radius lost, the signal likely needs a slightly more permissive peripyramidal field.
3. Error pattern: which donors are consistently wrong and what they share  
{summary_text}
{error_lines}

## Rationale
The best variation is biologically coherent because it asks whether the CA1 astrocyte-reactivity signal is concentrated specifically around CA1 pyramidal neurons rather than diffusely across all of CA1. That beat the nearby alternative because the winning radius better matched the scale at which astrocyte activation appears to covary with donor memory-decline burden. Since the current accepted panel already contains a region-wide astrocytic reactivity member, this niche-refined candidate is likely partly redundant but could still add value if the evaluator finds that neuron-adjacent gliosis carries cleaner or slightly more specific signal.

## Interpretation
The signal seems to represent CA1 peripyramidal gliosis: the population is CA1 astrocyte-lineage cells, the niche is immediate proximity to CA1 pyramidal neurons, the donor-level summary is the reactive fraction among those niche astrocytes, and the simplest observable pattern is more reactive astrocytes crowding the CA1 pyramidal field in donors with worse memory decline.

## Next
Run one more local sweep around the winning niche scale by testing a slightly tighter or slightly broader CA1 pyramidal-neighborhood threshold and, if needed, a denominator-stabilized variant such as a shrinkage-adjusted reactive fraction.
\"\"\"
"""
text = text.replace(old, new)
path.write_text(text)
PY
python -m py_compile /scratch/result.py
