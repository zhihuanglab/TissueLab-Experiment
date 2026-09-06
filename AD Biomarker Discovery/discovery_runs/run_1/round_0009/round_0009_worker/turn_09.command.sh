set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis.stats import leave_one_out_partial_correlation
print(inspect.signature(leave_one_out_partial_correlation))
print(inspect.getsource(leave_one_out_partial_correlation)[:4000])
PY
