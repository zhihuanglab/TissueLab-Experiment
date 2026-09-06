set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.append('/shared/lib')
from shared_analysis import stats
for fn in ['leave_one_out_summary','bootstrap_partial_correlation','leave_one_out_partial_correlation','residualize']:
    print(f'\n--- {fn} ---')
    print(inspect.getsource(getattr(stats, fn)))
PY
