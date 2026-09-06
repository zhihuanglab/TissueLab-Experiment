set -euo pipefail
python - <<'PY'
import sys, os
print('PYTHONPATH env=', os.environ.get('PYTHONPATH'))
print(sys.path)
PY
