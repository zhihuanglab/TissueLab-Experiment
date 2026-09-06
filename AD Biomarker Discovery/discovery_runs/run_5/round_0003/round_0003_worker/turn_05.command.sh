set -euo pipefail
echo 'shared dirs:'
find /shared -maxdepth 3 -type d | sed -n '1,60p'
echo 'python sys.path default:'
python - <<'PY'
import sys
print(sys.path)
PY
