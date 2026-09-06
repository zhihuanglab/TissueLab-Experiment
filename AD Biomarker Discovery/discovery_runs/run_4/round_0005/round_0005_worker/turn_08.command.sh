set -euo pipefail
find /data/autoresearch_runs -name results.json | head -n 5 | while read f; do
  echo "--- $f"
  python - <<PY
import json
f='$f'
try:
    d=json.load(open(f))
    print({k:d.get(k) for k in ['feature_name','best_variation','selection_score','loo_predictive_r','is_loo_gap','adjusted_score','penalty','partial_r'] if k in d})
    if 'ranked_variations' in d:
        print('ranked first', d['ranked_variations'][:2])
except Exception as e:
    print('ERR',e)
PY
done
