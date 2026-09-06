set -euo pipefail
python - <<'PY'
import json
with open('/scratch/results.json') as f:
    res=json.load(f)
best=[r for r in res['ranked_variations'] if r['variation_name']==res['best_variation']][0]
print(best)
PY
