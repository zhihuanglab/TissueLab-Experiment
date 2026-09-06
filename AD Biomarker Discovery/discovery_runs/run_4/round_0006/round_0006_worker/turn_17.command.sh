python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
txt=p.read_text()
txt=txt.replace('BEST_VARIATION = "reactive_pyramidal_desert_t10"', 'BEST_VARIATION = "reactive_pyramidal_desert_t20"')
txt=txt.replace('FEATURE_COLUMN = "ca1_reactive_pyramidal_desert_t10"', 'FEATURE_COLUMN = "ca1_reactive_pyramidal_desert_t20"')
txt=txt.replace(
'''    worked_why = (
        f"The winner, {best.name}, tracks the donor-level fraction of CA1 reactive astrocytes sitting in "
        f"pyramidal-neuron-poor microenvironments within a 35 µm CA1 niche. That is biologically aligned with a "
        f"reactive glial response to local neuron depletion or scar-like remodeling, so it matches the existing CA1 injury theme."
    )
''',
'''    worked_why = (
        f"The winner, {best.name}, produced the strongest in-sample residual association in this local sweep by tracking the donor-level "
        f"fraction of CA1 reactive astrocytes sitting in pyramidal-neuron-poor microenvironments within a 35 µm CA1 niche. "
        f"That biological framing is coherent with a reactive glial response to local neuron depletion or scar-like remodeling."
    )
''')
txt=txt.replace(
'''        failed_why = (
            f"The nearby alternative {runner_up.name} was weaker, suggesting that its threshold either labeled too many "
            f"reactive astrocytes as 'desert' cells or blurred the contrast between clearly neuron-depleted niches and more mixed neighborhoods."
        )
''',
'''        failed_why = (
            f"The nearby alternative {runner_up.name} was weaker, and more importantly the family as a whole generalized poorly "
            f"(winner LOO predictive r {best.loo_predictive_r:.4f}). That suggests the thresholded desert fraction is biologically sensible "
            f"but still too unstable or too redundant with existing CA1 injury features to yield reliable donor-level prediction on its own."
        )
''')
txt=txt.replace(
'''    rationale = (
        f"{best.name} won because the threshold of {VARIATIONS[best.name]['threshold']:.0%} keeps the biomarker focused on "
        f"reactive astrocytes occupying sharply neuron-depleted CA1 pockets, instead of merely reactive-rich but still neuron-mixed neighborhoods. "
        f"That makes it a coherent local refinement of the accepted panel's broader CA1 neuron-loss and reactive-astrocyte features."
    )
''',
'''    rationale = (
        f"{best.name} won because the threshold of {VARIATIONS[best.name]['threshold']:.0%} is less extreme than the 10% cutoff, "
        f"so it captures a broader band of neuron-poor reactive niches rather than only the most severe deserts. "
        f"That likely explains why it outperformed the nearby alternative locally, although the weak LOO signal means the family still looks unstable or redundant."
    )
''')
txt=txt.replace(
'''If relevant to the panel, this niche-focused reactive-astrocyte measure appears more specific than a global reactive burden feature and could plausibly add information beyond the current panel if it is not fully redundant with prior CA1 reactive-purity and hypertrophy members.
''',
'''Relative to the current panel, this niche-focused reactive-astrocyte measure is biologically more specific than a global reactive burden feature, but the near-zero LOO correlation makes it unlikely to add robust new information without another refinement.
''')
p.write_text(txt)
PY
python /scratch/result.py
