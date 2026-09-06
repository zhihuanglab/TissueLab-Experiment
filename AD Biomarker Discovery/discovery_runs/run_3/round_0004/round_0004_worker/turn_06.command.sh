set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.append('/shared/lib')
from shared_analysis.stats import leave_one_out_summary, partial_correlation
from shared_analysis.artifacts import build_results_payload, write_results_payload, validate_results_payload
print(inspect.getsource(partial_correlation))
print('---')
print(inspect.getsource(leave_one_out_summary))
print('---')
print(inspect.signature(build_results_payload))
print(inspect.signature(write_results_payload))
print(inspect.signature(validate_results_payload))
PY
