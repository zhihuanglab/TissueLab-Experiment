set -euo pipefail
export PYTHONPATH=/shared/lib
python - <<'PY'
import os, sys
print(os.environ.get('PYTHONPATH'))
print(sys.path[:5])
import shared_analysis
print(shared_analysis)
PY
