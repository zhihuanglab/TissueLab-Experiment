set -euo pipefail
echo '--- results.json summary ---'
python - <<'PY'
import json
with open('/scratch/results.json') as f:
    obj=json.load(f)
for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score','loo_unstable_count','loo_max_shift']:
    print(k, obj[k])
print('ranked len', len(obj['ranked_variations']))
print('best per_donor sample', obj['ranked_variations'][0]['per_donor_loo'][:2])
PY
echo '--- report.md ---'
sed -n '1,240p' /scratch/report.md
