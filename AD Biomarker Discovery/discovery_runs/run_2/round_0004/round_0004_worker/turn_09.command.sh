set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis.stats import partial_correlation, leave_one_out_summary, residualized_loo_predictive_correlation
from shared_analysis.artifacts import write_donor_feature_table, write_results_payload, build_results_payload
print('=== partial_correlation ===')
print(inspect.getsource(partial_correlation))
print('\n=== leave_one_out_summary ===')
print(inspect.getsource(leave_one_out_summary))
print('\n=== residualized_loo_predictive_correlation ===')
print(inspect.getsource(residualized_loo_predictive_correlation))
print('\n=== write_donor_feature_table ===')
print(inspect.getsource(write_donor_feature_table))
print('\n=== build_results_payload ===')
print(inspect.getsource(build_results_payload))
print('\n=== write_results_payload ===')
print(inspect.getsource(write_results_payload))
PY
