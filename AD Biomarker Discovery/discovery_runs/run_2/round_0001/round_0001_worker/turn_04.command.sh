set -euo pipefail
ls -la /shared
find /shared -maxdepth 2 -type d | sed -n '1,40p'
echo 'PYTHONPATH='$PYTHONPATH
python - <<'PY'
import sys, os
print(sys.path[:10])
PY
