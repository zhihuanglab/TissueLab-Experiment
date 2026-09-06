set -euo pipefail
python - <<'PY'
import importlib.util, json
from pathlib import Path
spec = importlib.util.spec_from_file_location('worker_result', '/scratch/result.py')
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
print('canonical variation:', mod.CANONICAL_REPLAY_VARIATION)
print('sample score:', mod.compute_donor_score(donor_id='H19.33.004', data_root='/data'))
print('results keys:', sorted(json.loads(Path('/scratch/results.json').read_text()).keys())[:20])
print('report exists:', Path('/scratch/report.md').exists())
print('donor table exists:', Path('/scratch/donor_feature_table.csv').exists())
PY
