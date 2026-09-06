set -euo pipefail
find /data/autoresearch_runs -name results.json | head -5 | while read f; do echo '---' $f; python - <<PY
import json,sys
with open("$f") as fh:
    j=json.load(fh)
for k in ['best_variation','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score','feature_column','ranked_variations']:
    if k in j:
        print(k, j[k] if k!='ranked_variations' else j[k][:2])
PY
done
