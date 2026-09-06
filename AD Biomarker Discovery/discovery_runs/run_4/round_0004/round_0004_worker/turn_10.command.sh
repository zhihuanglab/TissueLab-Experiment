set -euo pipefail
python - <<'PY'
from pathlib import Path
for p in [Path('/shared/templates/worker_analysis_template.py'), Path('/shared/templates/worker_embedding_mechanistic_template.py')]:
    print('---', p, '---')
    txt=p.read_text()
    for needle in ['Adjusted Score','selection_score','is_loo_gap']:
        if needle in txt:
            print('contains', needle)
    lines=txt.splitlines()
    for i,l in enumerate(lines):
        if 'selection_score' in l or 'adjusted_score' in l or 'loo_predictive' in l:
            print(i+1, l)
PY
