set -euo pipefail
python - <<'PY'
from pathlib import Path
import json, glob
paths=sorted(glob.glob('/data/autoresearch_runs/**/*/results.json', recursive=True))
print('n results', len(paths))
for p in paths[-3:]:
    print('---', p)
    with open(p) as f:
        d=json.load(f)
    print({k:d.get(k) for k in ['best_variation','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score','feature_column']})
    if 'ranked_variations' in d:
        print('ranked first', d['ranked_variations'][:2])
PY
