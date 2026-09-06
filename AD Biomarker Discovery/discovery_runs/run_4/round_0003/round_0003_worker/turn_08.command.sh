set -e
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis.stats import partial_correlation, leave_one_out_summary, residualized_loo_predictive_correlation
print('partial_correlation\n', inspect.getsource(partial_correlation))
print('\nleave_one_out_summary head\n', '\n'.join(inspect.getsource(leave_one_out_summary).splitlines()[:220]))
print('\nresidualized_loo_predictive_correlation\n', inspect.getsource(residualized_loo_predictive_correlation))
PY
