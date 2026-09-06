set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
from shared_analysis.stats import partial_correlation, leave_one_out_summary
print(inspect.getsource(partial_correlation))
print('---')
print(inspect.getsource(leave_one_out_summary))
PY
