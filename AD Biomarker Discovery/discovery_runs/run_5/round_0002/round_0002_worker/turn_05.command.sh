set -euo pipefail
python - <<'PY'
import sys, inspect, textwrap
sys.path.insert(0, '/shared/lib')
from shared_analysis import partial_correlation, residualized_loo_predictive_correlation, write_donor_feature_table, write_results_payload
print('partial_correlation:\n', textwrap.dedent(inspect.getsource(partial_correlation)))
print('\nresidualized_loo_predictive_correlation:\n', textwrap.dedent(inspect.getsource(residualized_loo_predictive_correlation)))
print('\nwrite_donor_feature_table:\n', textwrap.dedent(inspect.getsource(write_donor_feature_table)))
print('\nwrite_results_payload:\n', textwrap.dedent(inspect.getsource(write_results_payload)))
PY
