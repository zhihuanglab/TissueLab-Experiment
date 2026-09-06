set -euo pipefail
PYTHONPATH=/shared/lib python3 - <<'PY'
import os, sys
print(os.environ.get('PYTHONPATH'))
print(sys.path[:6])
import shared_analysis
print(shared_analysis.__file__)
PY
