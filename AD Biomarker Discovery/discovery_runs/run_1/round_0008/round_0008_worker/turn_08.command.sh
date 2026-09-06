set -euo pipefail
find /data/autoresearch_runs -name results.json | head -n 3
python - <<'PY'
import json, glob
paths=glob.glob('/data/autoresearch_runs/**/results.json', recursive=True)
print('n', len(paths))
for p in paths[:3]:
    print('---', p)
    with open(p) as f:
        obj=json.load(f)
    for k in ['best_variation','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score','feature_column']:
        if k in obj: print(k, obj[k])
    if 'ranked_variations' in obj:
        print('first var', obj['ranked_variations'][0])
PY
