set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
from shared_analysis.stats import partial_correlation, residualized_loo_predictive_correlation, leave_one_out_summary
print('partial_correlation', inspect.signature(partial_correlation))
print(inspect.getsource(partial_correlation)[:2200])
print('resid loo', inspect.signature(residualized_loo_predictive_correlation))
print(inspect.getsource(residualized_loo_predictive_correlation)[:2600])
print('loo summary', inspect.signature(leave_one_out_summary))
print(inspect.getsource(leave_one_out_summary)[:2200])
PY
