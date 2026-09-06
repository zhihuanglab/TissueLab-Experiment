set -euo pipefail
echo '--- report.md ---'
sed -n '1,240p' /scratch/report.md
echo '--- results.json head ---'
python - <<'PY'
import json
obj=json.load(open('/scratch/results.json'))
for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score']:
    print(k, obj.get(k))
print('ranked len', len(obj.get('ranked_variations',[])))
PY
