set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis.stats import partial_correlation, residualized_loo_predictive_correlation
from shared_analysis.artifacts import write_donor_feature_table, write_results_payload, build_results_payload
print('partial_correlation:', inspect.signature(partial_correlation))
print(inspect.getsource(partial_correlation)[:3500])
print('--- residualized_loo_predictive_correlation:', inspect.signature(residualized_loo_predictive_correlation))
print(inspect.getsource(residualized_loo_predictive_correlation)[:4500])
print('--- write_donor_feature_table:', inspect.signature(write_donor_feature_table))
print(inspect.getsource(write_donor_feature_table)[:2500])
print('--- build_results_payload:', inspect.signature(build_results_payload))
print(inspect.getsource(build_results_payload)[:3000])
PY
