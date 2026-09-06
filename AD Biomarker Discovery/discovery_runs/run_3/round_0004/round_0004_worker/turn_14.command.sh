set -euo pipefail
sed -n '1,240p' /scratch/report.md
echo
echo '--- results.json head ---'
python - <<'PY'
import json
from pathlib import Path
d=json.loads(Path('/scratch/results.json').read_text())
print(json.dumps({k:d[k] for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','gap_penalty','adjusted_score']}, indent=2))
print('ranked names:', [x['variation_name'] for x in d['ranked_variations']])
PY
