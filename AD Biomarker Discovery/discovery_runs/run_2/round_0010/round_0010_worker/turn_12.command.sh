set -e
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis import bootstrap_partial_correlation, leave_one_out_summary, build_results_payload, write_donor_feature_table
from shared_analysis.stats import residualized_loo_predictive_correlation, bootstrap_partial_correlation_stability
print('bootstrap_partial_correlation:', inspect.signature(bootstrap_partial_correlation))
print(inspect.getsource(bootstrap_partial_correlation)[:1500])
print('leave_one_out_summary:', inspect.signature(leave_one_out_summary))
print(inspect.getsource(leave_one_out_summary)[:1500])
print('residualized_loo_predictive_correlation:', inspect.signature(residualized_loo_predictive_correlation))
print(inspect.getsource(residualized_loo_predictive_correlation)[:2000])
print('write_donor_feature_table:', inspect.signature(write_donor_feature_table))
print(inspect.getsource(write_donor_feature_table)[:1500])
PY
