set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis.stats import residualized_loo_predictive_correlation, resolve_covariate_names, leave_one_out_partial_correlation
print('residualized_loo_predictive_correlation sig:', inspect.signature(residualized_loo_predictive_correlation))
print(inspect.getsource(residualized_loo_predictive_correlation))
print('\nleave_one_out_partial_correlation sig:', inspect.signature(leave_one_out_partial_correlation))
print(inspect.getsource(leave_one_out_partial_correlation)[:4000])
PY
