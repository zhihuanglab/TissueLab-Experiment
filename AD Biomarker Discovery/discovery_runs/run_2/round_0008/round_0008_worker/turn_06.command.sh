set -euo pipefail
echo '--- round 7 result.py ---'
sed -n '1,280p' /data/autoresearch_runs/run_20260420_215745_90b105/round_0007/round_0007_worker/sandbox/result.py
echo
echo '--- round 7 results.json excerpt ---'
python - <<'PY'
import json
path='/data/autoresearch_runs/run_20260420_215745_90b105/round_0007/round_0007_worker/sandbox/results.json'
obj=json.load(open(path))
for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score']:
    print(k, obj.get(k))
print('ranked variations:')
for row in obj.get('ranked_variations',[]):
    print(row)
PY
