set -euo pipefail
python - <<'PY'
import sys, os
print(sys.executable)
print('PYTHONPATH env', os.environ.get('PYTHONPATH'))
print(sys.path[:5])
PY
