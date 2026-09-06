set -euo pipefail
python - <<'PY'
import os, json, glob
paths=glob.glob('/data/autoresearch_runs/**/*results.json', recursive=True)[:10]
print('found', len(paths))
for p in paths[:5]:
    try:
        obj=json.load(open(p))
    except Exception:
        continue
    if 'adjusted_score' in obj or 'selection_score' in obj:
        print(p)
        print({k:obj.get(k) for k in ['feature_name','selection_score','loo_predictive_r','is_loo_gap','adjusted_score','best_variation']})
        break
PY
