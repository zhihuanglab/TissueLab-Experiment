set -euo pipefail
echo 'searching prior artifacts for "Adjusted Score" or "penalty="'
grep -R "Adjusted Score\|penalty=" -n /data/autoresearch_runs 2>/dev/null | head -40 || true
echo 'searching previous results.json examples'
find /data/autoresearch_runs -name results.json | head -5 | while read f; do echo "--- $f"; python - <<PY
import json; import sys
p=sys.argv[1]
obj=json.load(open(p))
for k in ('best_variation','ranked_variations','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score','feature_column','feature_name'):
    if k in obj: print(k, obj[k])
PY "$f"; done
