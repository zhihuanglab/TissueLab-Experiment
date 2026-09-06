set -euo pipefail
echo 'Shared tree'
find /shared -maxdepth 3 -type f | sed -n '1,120p'
echo 'Python version and sys.path'
python - <<'PY'
import sys
print(sys.version)
print(sys.path)
PY
