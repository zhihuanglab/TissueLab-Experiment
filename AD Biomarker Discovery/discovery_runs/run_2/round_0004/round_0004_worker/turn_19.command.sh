set -euo pipefail
python - <<'PY'
import json, pandas as pd
res=json.loads(open('/scratch/results.json').read())
print(pd.DataFrame(res['ranked_variations'])[['variation_name','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score','feature_column']].to_string(index=False))
PY
