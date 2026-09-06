set -euo pipefail
echo '--- canonical target in result.py ---'
grep -n "BEST_VARIATION_NAME\|BEST_VARIATION_RADIUS_PX\|BEST_FEATURE_COLUMN" /scratch/result.py
echo '--- results.json head ---'
python - <<'PY'
import json
with open('/scratch/results.json') as f:
    obj=json.load(f)
print(json.dumps({
    'best_variation': obj.get('best_variation'),
    'feature_column': obj.get('feature_column'),
    'n_total': obj.get('n_total'),
    'n_analyzable': obj.get('n_analyzable'),
    'partial_r': obj.get('partial_r'),
    'selection_score': obj.get('selection_score'),
    'loo_predictive_r': obj.get('loo_predictive_r'),
    'ranked_variations': obj.get('ranked_variations'),
}, indent=2)[:6000])
PY
echo '--- report.md ---'
sed -n '1,220p' /scratch/report.md
