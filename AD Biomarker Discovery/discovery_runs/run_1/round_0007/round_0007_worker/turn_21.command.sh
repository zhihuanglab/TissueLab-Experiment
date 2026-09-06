set -euo pipefail
python /scratch/result.py > /scratch/result_stdout.txt
tail -n 160 /scratch/result_stdout.txt
echo '--- results.json head ---'
python - <<'PY'
import json
with open('/scratch/results.json') as f:
    obj=json.load(f)
print(json.dumps({k: obj[k] for k in ['best_variation','feature_name','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','gap_penalty','adjusted_score']}, indent=2))
print('ranked count', len(obj['ranked_variations']))
PY
echo '--- report.md ---'
sed -n '1,220p' /scratch/report.md
