set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0,'/shared/lib')
from shared_analysis.artifacts import _COVARIATE_ALIASES
print(_COVARIATE_ALIASES)
PY
