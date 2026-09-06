set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
python - <<'PY'
import inspect
from shared_analysis import build_cell_table, load_training_cohort
print('load_training_cohort:\n', inspect.getsource(load_training_cohort))
print('\n--- build_cell_table signature ---')
print(inspect.signature(build_cell_table))
print('\n--- build_cell_table source head ---')
src = inspect.getsource(build_cell_table)
print(src[:8000])
PY
