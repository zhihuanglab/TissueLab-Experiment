set -euo pipefail
python /scratch/result.py > /scratch/run_stdout.txt
cat /scratch/run_stdout.txt
echo '--- results.json head ---'
python - <<'PY'
import json
from pathlib import Path
p=Path('/scratch/results.json')
print('exists', p.exists(), 'size', p.stat().st_size if p.exists() else None)
if p.exists():
    obj=json.loads(p.read_text())
    print({k: obj.get(k) for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score']})
    print('ranked count', len(obj.get('ranked_variations', [])))
PY
echo '--- report.md head ---'
sed -n '1,220p' /scratch/report.md
