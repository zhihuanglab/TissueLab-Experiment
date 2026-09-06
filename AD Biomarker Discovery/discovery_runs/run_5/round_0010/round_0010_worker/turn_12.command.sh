set -euo pipefail
echo '--- results.json ---'
python - <<'PY'
import json
data=json.load(open('/scratch/results.json'))
print(json.dumps({k:data[k] for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score']}, indent=2))
print('ranked_variations len', len(data['ranked_variations']))
PY
echo '--- report.md ---'
sed -n '1,240p' /scratch/report.md
