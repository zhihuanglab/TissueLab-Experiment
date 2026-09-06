set -euo pipefail
find /shared -maxdepth 3 -type d | sed -n '1,120p'
echo 'PYTHONPATH='$PYTHONPATH
python - <<'PY'
import sys, os
print(sys.path)
PY
