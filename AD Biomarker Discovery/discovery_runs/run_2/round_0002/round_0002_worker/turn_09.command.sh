set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis.stats import partial_correlation, leave_one_out_summary, residualized_loo_predictive_correlation
print(inspect.getsource(partial_correlation))
print('\n--- residualized_loo_predictive_correlation ---')
print(inspect.getsource(residualized_loo_predictive_correlation))
print('\n--- leave_one_out_summary ---')
print(inspect.getsource(leave_one_out_summary))
PY
