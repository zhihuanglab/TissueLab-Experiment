set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.append('/shared/lib')
from shared_analysis.stats import residualized_loo_predictive_correlation, leave_one_out_partial_correlation
print(inspect.getsource(residualized_loo_predictive_correlation))
print('---')
print(inspect.getsource(leave_one_out_partial_correlation))
PY
