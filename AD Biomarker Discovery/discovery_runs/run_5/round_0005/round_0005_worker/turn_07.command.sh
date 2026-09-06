set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.append('/shared/lib')
from shared_analysis.artifacts import build_results_payload, write_results_payload, coerce_results_payload, validate_results_payload
print(inspect.getsource(build_results_payload))
print(inspect.getsource(validate_results_payload))
PY
