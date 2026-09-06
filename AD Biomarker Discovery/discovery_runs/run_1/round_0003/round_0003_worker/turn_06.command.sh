set -e
PYTHONPATH=/shared/lib python - <<'PY'
import sys, os
print(os.environ.get('PYTHONPATH'))
print(sys.path[:10])
import importlib.util
print(importlib.util.find_spec('shared_analysis'))
PY
