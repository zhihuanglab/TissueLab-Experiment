set -e
python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
text=p.read_text()
text=text.replace('            "cognitive_status": str(getattr(row, "cognitive_status")),\n','')
p.write_text(text)
PY
python /scratch/result.py
echo
echo '--- results.json ---'
python - <<'PY'
import json
obj=json.load(open('/scratch/results.json'))
print(json.dumps({k: obj[k] for k in ['best_variation','feature_column','selection_score','partial_r','loo_predictive_r','n_analyzable','n_total']}, indent=2))
print('ranked_variations:', [x['name'] for x in obj['ranked_variations']])
PY
echo
echo '--- report.md head ---'
sed -n '1,220p' /scratch/report.md
