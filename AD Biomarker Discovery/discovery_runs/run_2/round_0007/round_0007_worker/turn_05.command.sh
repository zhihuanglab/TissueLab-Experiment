set -euo pipefail
echo 'PYTHONPATH='${PYTHONPATH:-}
find /shared -maxdepth 3 -type d | sed -n '1,80p'
find /shared -maxdepth 4 -type f | grep -E 'shared_analysis|template|runtime_quickstart' | sed -n '1,80p'
python - <<'PY'
import sys, os
print(sys.path)
PY
