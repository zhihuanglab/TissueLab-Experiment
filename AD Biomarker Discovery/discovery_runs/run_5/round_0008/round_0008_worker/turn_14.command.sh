set -euo pipefail
echo '--- updated canonical line ---'
grep -n 'CANONICAL_VARIATION' /scratch/result.py | head
echo
echo '--- report.md ---'
sed -n '1,220p' /scratch/report.md
echo
echo '--- results.json head ---'
python - <<'PY'
import json
with open('/scratch/results.json') as f: obj=json.load(f)
print(json.dumps({k:obj[k] for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score']}, indent=2))
PY
