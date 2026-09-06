set -euo pipefail
python - <<'PY'
import sys, os
print(sys.version)
print('PYTHONPATH', os.environ.get('PYTHONPATH'))
print('sys.path head', sys.path[:10])
PY
echo 'shared lib listing'
find /shared -maxdepth 3 -type f | sed -n '1,80p'
