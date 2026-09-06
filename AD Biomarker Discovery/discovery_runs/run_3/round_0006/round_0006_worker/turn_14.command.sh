set -euo pipefail
python - <<'PY'
import sys, inspect, textwrap
sys.path.insert(0,'/shared/lib')
from shared_analysis.artifacts import resolve_covariate_names, normalize_covariate_names
print(textwrap.dedent(inspect.getsource(resolve_covariate_names)))
print(textwrap.dedent(inspect.getsource(normalize_covariate_names)))
PY
