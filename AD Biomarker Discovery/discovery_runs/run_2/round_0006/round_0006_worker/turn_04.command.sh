set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
python - <<'PY'
from pathlib import Path
import inspect
from shared_analysis import build_cell_table, load_training_cohort
print('build_cell_table signature:', inspect.signature(build_cell_table))
print('\nsource build_cell_table:\n')
src=inspect.getsource(build_cell_table)
print(src[:7000])
print('\nload_training_cohort sig:', inspect.signature(load_training_cohort))
PY
echo '--- previous round 5 result.py ---'
sed -n '1,260p' /data/autoresearch_runs/run_20260420_215745_90b105/round_0005/round_0005_worker/result.py || true
echo '--- previous round 3 result.py ---'
sed -n '1,260p' /data/autoresearch_runs/run_20260420_215745_90b105/round_0003/round_0003_worker/result.py || true
