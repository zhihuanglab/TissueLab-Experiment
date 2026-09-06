set -euo pipefail
python - <<'PY'
import json, glob
paths=sorted(glob.glob('/data/autoresearch_runs/**/round_0006/**/results.json', recursive=True))
for p in paths[:2]:
    print('---', p)
    with open(p) as f:
        obj=json.load(f)
    for k in ['best_variation','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score','feature_column','feature_name']:
        print(k, obj.get(k))
    print('ranked vars head:')
    for rv in obj.get('ranked_variations',[])[:3]:
        print(rv)
PY
