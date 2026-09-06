set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0,'/shared/lib')
from shared_analysis import DEFAULT_CONFOUNDS
print(DEFAULT_CONFOUNDS)
PY
