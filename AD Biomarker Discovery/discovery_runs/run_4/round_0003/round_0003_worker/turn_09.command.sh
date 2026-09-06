set -e
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis.artifacts import write_donor_feature_table, write_results_payload, build_results_payload, validate_results_payload
print('write_donor_feature_table\n', inspect.getsource(write_donor_feature_table))
print('\nbuild_results_payload head\n', '\n'.join(inspect.getsource(build_results_payload).splitlines()[:220]))
print('\nvalidate_results_payload head\n', '\n'.join(inspect.getsource(validate_results_payload).splitlines()[:180]))
PY
