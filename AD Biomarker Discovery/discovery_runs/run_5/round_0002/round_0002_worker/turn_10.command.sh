set -euo pipefail
python - <<'PY'
import sys, inspect, textwrap
sys.path.insert(0, '/shared/lib')
from shared_analysis.artifacts import build_results_payload, validate_results_payload
print(textwrap.dedent(inspect.getsource(build_results_payload)))
print('--- validate ---')
print(textwrap.dedent(inspect.getsource(validate_results_payload)))
PY
