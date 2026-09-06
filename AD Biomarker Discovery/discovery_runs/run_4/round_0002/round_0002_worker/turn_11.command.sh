set -e
echo '--- result.py head ---'
sed -n '1,80p' /scratch/result.py
echo
echo '--- results.json ---'
python - <<'PY'
import json
p='/scratch/results.json'
data=json.load(open(p))
print(json.dumps({k:data[k] for k in ['best_variation','feature_column','selection_score','partial_r','loo_predictive_r','is_loo_gap','penalty','adjusted_score']}, indent=2))
print('ranked count', len(data['ranked_variations']))
PY
echo
echo '--- report.md ---'
sed -n '1,260p' /scratch/report.md
