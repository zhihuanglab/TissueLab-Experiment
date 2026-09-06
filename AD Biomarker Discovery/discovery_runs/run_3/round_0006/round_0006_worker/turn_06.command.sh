set -euo pipefail
PYTHONPATH=/shared/lib python - <<'PY'
import sys
print(sys.path[:5])
import shared_analysis
print(shared_analysis)
PY
