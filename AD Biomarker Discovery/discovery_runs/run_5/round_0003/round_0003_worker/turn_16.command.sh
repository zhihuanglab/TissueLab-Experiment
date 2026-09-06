python - <<'PY'
import json
from pathlib import Path
data=json.loads(Path('/scratch/results.json').read_text())
print(json.dumps({k:data[k] for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score']}, indent=2))
PY
