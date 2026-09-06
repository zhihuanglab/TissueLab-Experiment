set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
import shared_analysis
print(shared_analysis.DEFAULT_CONFOUNDS)
print(inspect.getsource(shared_analysis))
PY
