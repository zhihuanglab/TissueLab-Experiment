set -euo pipefail
python3 - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
old = """    analyzable = int(best["n_analyzable"])
    n_total = int(best["n_total"])

    return f\"\"\"## Summary
"""
new = """    analyzable = int(best["n_analyzable"])
    n_total = int(best["n_total"])

    if loser_rows:
        if float(best["radius_px"]) < float(loser["radius_px"]):
            rationale_beat_text = (
                "It beat the nearby alternative because the winning radius stayed tighter around the "
                "perineuronal zone and avoided diluting the signal with more diffuse background CA1 "
                "reactive astrocyte burden."
            )
        else:
            rationale_beat_text = (
                "It beat the nearby alternative because the winning radius retained enough short-range "
                "context to capture neuron-associated astrocyte exposure without becoming so strict "
                "that genuine perineuronal coverage events were missed."
            )
    else:
        rationale_beat_text = "It beat the nearby alternative because no alternate local variation outperformed it."

    return f\"\"\"## Summary
"""
if old not in text:
    raise SystemExit('anchor not found')
text = text.replace(old, new)

old2 = """The best variation is biologically coherent because it conditions first on the vulnerable CA1 pyramidal neuron population and then measures local reactive astrocyte proximity around those neurons, rather than counting reactive astrocytes globally. That makes the scalar interpretable as peripyramidal reactive astrocyte coverage. It beat the nearby alternative because the winning radius retained enough local context to capture reactive-neighbor exposure while avoiding a too-stringent definition that likely undercounts genuine neuron-associated astrocytosis. Relative to the accepted panel, this candidate is plausibly somewhat additive if it captures conditional local exposure beyond overall CA1 pyramidal loss, reactive astrocyte burden, and the broader mixed-cell niche fraction already in panel review.
"""
new2 = """The best variation is biologically coherent because it conditions first on the vulnerable CA1 pyramidal neuron population and then measures local reactive astrocyte proximity around those neurons, rather than counting reactive astrocytes globally. That makes the scalar interpretable as peripyramidal reactive astrocyte coverage. {rationale_beat_text} Relative to the accepted panel, this candidate is plausibly somewhat additive if it captures conditional local exposure beyond overall CA1 pyramidal loss, reactive astrocyte burden, and the broader mixed-cell niche fraction already in panel review.
"""
if old2 not in text:
    raise SystemExit('rationale paragraph not found')
text = text.replace(old2, new2)
path.write_text(text, encoding='utf-8')
print('patched rationale text')
PY
/usr/local/bin/python3 /scratch/result.py >/scratch/run_stdout.txt
grep -n "## Rationale" -A4 /scratch/report.md
