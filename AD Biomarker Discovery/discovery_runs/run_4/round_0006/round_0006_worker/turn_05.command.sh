set -euo pipefail
echo "PYTHONPATH=${PYTHONPATH:-}"
ls -la /shared
find /shared -maxdepth 3 -type d | sed -n '1,60p'
python - <<'PY'
import sys, os
print(sys.path)
PY
