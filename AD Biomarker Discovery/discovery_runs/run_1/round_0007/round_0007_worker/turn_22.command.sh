python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()

text = text.replace('DEFAULT_CANONICAL_VARIATION = "candidate_variant_a"', 'DEFAULT_CANONICAL_VARIATION = "candidate_variant_b"')

old = '''    findings_failed = []
    for entry in ranked[1:]:
        findings_failed.append(
            f"{entry['name']} was weaker than the winner (selection {entry['selection_score']:.4f} vs {best['selection_score']:.4f}), "
            f"consistent with its broader small-cell cutoff diluting the niche-specific small-cell tail."
        )
    failed_text = " ".join(findings_failed) if findings_failed else "No nearby alternatives were tested."

    additivity_text = (
        "Because the accepted panel already contains CA1 pyramidal abundance, reactive astrocyte lineage burden, and a proximal median-area term, "
        "this feature is most plausible as a thresholded tail version of the same niche morphology rather than a wholly orthogonal biology."
    )

    return f"""## Summary
'''
new = '''    findings_failed = []
    for entry in ranked[1:]:
        if (entry.get("small_quantile") or 0.0) < (best.get("small_quantile") or 0.0):
            why = "its stricter tail cutoff may have been too sparse, missing moderately shrunken proximal pyramidal neurons"
        elif (entry.get("small_quantile") or 0.0) > (best.get("small_quantile") or 0.0):
            why = "its broader cutoff may have diluted the niche-specific small-cell burden with weakly shifted cells"
        else:
            why = "the nearby parameter change did not sharpen the niche definition"
        findings_failed.append(
            f"{entry['name']} was weaker than the winner (selection {entry['selection_score']:.4f} vs {best['selection_score']:.4f}), consistent with {why}."
        )
    failed_text = " ".join(findings_failed) if findings_failed else "No nearby alternatives were tested."

    if any((entry.get("small_quantile") or 0.0) < (best.get("small_quantile") or 0.0) for entry in ranked[1:]):
        quantile_rationale = "Using the 33rd percentile beat the stricter 25th percentile cutoff, which suggests the signal is not confined to the most extreme atrophic tail and is better captured by a somewhat broader burden of shrunken proximal cells."
    elif any((entry.get("small_quantile") or 0.0) > (best.get("small_quantile") or 0.0) for entry in ranked[1:]):
        quantile_rationale = "Using the winning cutoff outperformed broader alternatives, implying the cleanest signal comes from the more extreme small-cell tail rather than a mild left-shift of the whole proximal distribution."
    else:
        quantile_rationale = "The winning cutoff gave the strongest local single-feature signal in this family-local sweep."

    additivity_text = (
        "Because the accepted panel already contains CA1 pyramidal abundance, reactive astrocyte lineage burden, and a proximal median-area term, "
        "this feature is most plausible as a thresholded tail version of the same niche morphology rather than a wholly orthogonal biology."
    )

    return f"""## Summary
'''
if old not in text:
    raise SystemExit('first target block not found')
text = text.replace(old, new)

old2 = '''## Rationale
The best variation keeps the previously productive CA1 pyramidal–reactive astrocyte niche but converts the continuous area shift into a donor-normalized tail burden. Using the 25th percentile beat the broader 33rd percentile cutoff, which implies the strongest signal comes from the more extreme small-cell tail rather than a mild left-shift of the whole proximal distribution. {additivity_text}
'''
new2 = '''## Rationale
The best variation keeps the previously productive CA1 pyramidal–reactive astrocyte niche but converts the continuous area shift into a donor-normalized tail burden. {quantile_rationale} {additivity_text}
'''
if old2 not in text:
    raise SystemExit('second target block not found')
text = text.replace(old2, new2)

path.write_text(text)
PY
python -m py_compile /scratch/result.py
python - <<'PY'
import json
from pathlib import Path
res = json.loads(Path('/scratch/results.json').read_text())
best = res['best_variation']
ranked = res['ranked_variations']
best_entry = next(v for v in ranked if v['name']==best)
other = [v for v in ranked if v['name']!=best_entry['name']]
if other and other[0]['small_quantile'] < best_entry['small_quantile']:
    failed = f"{other[0]['name']} was weaker than the winner (selection {other[0]['selection_score']:.4f} vs {best_entry['selection_score']:.4f}), consistent with its stricter tail cutoff being too sparse and missing moderately shrunken proximal pyramidal neurons."
    rationale = "Using the 33rd percentile beat the stricter 25th percentile cutoff, which suggests the signal is not confined to the most extreme atrophic tail and is better captured by a somewhat broader burden of shrunken proximal cells."
else:
    failed = f"{other[0]['name']} was weaker than the winner (selection {other[0]['selection_score']:.4f} vs {best_entry['selection_score']:.4f})."
    rationale = "The winning cutoff gave the strongest local single-feature signal in this family-local sweep."
report = Path('/scratch/report.md').read_text()
report = report.replace(
"2. What failed and why (specific to the chosen hypothesis and what went wrong): candidate_variant_a was weaker than the winner (selection 0.2568 vs 0.3287), consistent with its broader small-cell cutoff diluting the niche-specific small-cell tail.",
f"2. What failed and why (specific to the chosen hypothesis and what went wrong): {failed}"
)
report = report.replace(
"The best variation keeps the previously productive CA1 pyramidal–reactive astrocyte niche but converts the continuous area shift into a donor-normalized tail burden. Using the 25th percentile beat the broader 33rd percentile cutoff, which implies the strongest signal comes from the more extreme small-cell tail rather than a mild left-shift of the whole proximal distribution. Because the accepted panel already contains CA1 pyramidal abundance, reactive astrocyte lineage burden, and a proximal median-area term, this feature is most plausible as a thresholded tail version of the same niche morphology rather than a wholly orthogonal biology.",
f"The best variation keeps the previously productive CA1 pyramidal–reactive astrocyte niche but converts the continuous area shift into a donor-normalized tail burden. {rationale} Because the accepted panel already contains CA1 pyramidal abundance, reactive astrocyte lineage burden, and a proximal median-area term, this feature is most plausible as a thresholded tail version of the same niche morphology rather than a wholly orthogonal biology."
)
Path('/scratch/report.md').write_text(report)
PY
sed -n '1,220p' /scratch/report.md
