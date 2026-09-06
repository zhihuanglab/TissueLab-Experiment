set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis import partial_correlation, residualized_loo_predictive_correlation, leave_one_out_summary
import shared_analysis.stats as s
print('partial_correlation sig:', inspect.signature(s.partial_correlation))
print(inspect.getsource(s.partial_correlation))
print('\nresidualized_loo_predictive_correlation sig:', inspect.signature(s.residualized_loo_predictive_correlation))
print(inspect.getsource(s.residualized_loo_predictive_correlation))
print('\nleave_one_out_summary sig:', inspect.signature(s.leave_one_out_summary))
print(inspect.getsource(s.leave_one_out_summary))
PY
