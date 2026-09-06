set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.append('/shared/lib')
import shared_analysis
print(shared_analysis.DEFAULT_CONFOUNDS)
PY
