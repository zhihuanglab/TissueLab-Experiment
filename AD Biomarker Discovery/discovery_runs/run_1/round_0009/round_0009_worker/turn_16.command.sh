set -euo pipefail
echo '--- results.json ---'
python - <<'PY'
import json
from pathlib import Path
p=Path('/scratch/results.json')
obj=json.loads(p.read_text())
for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score']:
    print(k, obj.get(k))
print('ranked len', len(obj.get('ranked_variations',[])))
print('top ranked', obj['ranked_variations'][0])
PY
echo '--- report.md ---'
sed -n '1,240p' /scratch/report.md
echo '--- canonical markers ---'
grep -n "CANONICAL_VARIATION\\|FEATURE_NAME\\|FEATURE_COLUMN" /scratch/result.py | sed -n '1,20p'
