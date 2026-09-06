set -euo pipefail
echo '--- results.json ---'
python - <<'PY'
import json, pathlib
p=pathlib.Path('/scratch/results.json')
obj=json.loads(p.read_text())
for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score']:
    print(k, obj.get(k))
print('ranked variations n=', len(obj.get('ranked_variations', [])))
PY
echo '--- report.md ---'
sed -n '1,220p' /scratch/report.md
echo '--- result_state.json ---'
cat /scratch/result_state.json
