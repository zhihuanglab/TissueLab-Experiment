set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.append('/shared/lib')
from shared_analysis.stats import partial_correlation, leave_one_out_summary, residualized_loo_predictive_correlation
from shared_analysis.artifacts import write_donor_feature_table, build_results_payload, write_results_payload
print('partial_correlation:\n', inspect.getsource(partial_correlation))
print('residualized_loo_predictive_correlation:\n', inspect.getsource(residualized_loo_predictive_correlation))
print('leave_one_out_summary:\n', inspect.getsource(leave_one_out_summary))
print('write_donor_feature_table:\n', inspect.getsource(write_donor_feature_table))
PY
