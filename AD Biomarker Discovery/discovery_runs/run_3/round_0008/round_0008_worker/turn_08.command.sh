set -euo pipefail
python - <<'PY'
import sys, inspect, json
sys.path.append('/shared/lib')
import shared_analysis
from shared_analysis import stats, artifacts
print('shared_analysis exports:', [a for a in dir(shared_analysis) if not a.startswith('_')][:50])
print('\npartial_correlation:\n', inspect.getsource(stats.partial_correlation))
print('\nresidualized_loo_predictive_correlation:\n', inspect.getsource(stats.residualized_loo_predictive_correlation))
print('\nwrite_donor_feature_table:\n', inspect.getsource(artifacts.write_donor_feature_table))
print('\nbuild_results_payload:\n', inspect.getsource(artifacts.build_results_payload))
PY
