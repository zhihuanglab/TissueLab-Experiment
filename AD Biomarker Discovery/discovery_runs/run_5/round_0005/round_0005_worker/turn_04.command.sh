set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.append('/shared/lib')
from shared_analysis import build_cell_table, load_training_cohort
print('loaded shared_analysis ok')
print(inspect.getsource(build_cell_table))
PY
