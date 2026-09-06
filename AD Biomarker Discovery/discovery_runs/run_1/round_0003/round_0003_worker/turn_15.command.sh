set -e
ls -l /scratch/donor_feature_table.csv /scratch/results.json /scratch/result.py /scratch/run_stdout.txt
python - <<'PY'
import json
from pathlib import Path
res = json.loads(Path('/scratch/results.json').read_text())
print(json.dumps({k:res[k] for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score']}, indent=2))
print('ranked len', len(res['ranked_variations']))
PY
