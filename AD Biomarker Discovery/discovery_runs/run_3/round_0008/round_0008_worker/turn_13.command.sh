set -euo pipefail
python - <<'PY'
import json, glob, os
for p in glob.glob('/data/autoresearch_runs/**/results.json', recursive=True):
    try:
        obj=json.load(open(p))
    except Exception:
        continue
    keys=set(obj)
    if 'selection_score' in keys:
        print('path', p)
        print({k:obj.get(k) for k in ['feature_name','selection_score','loo_predictive_r','is_loo_gap','adjusted_score','partial_r','n_analyzable','n_total']})
        break
PY
