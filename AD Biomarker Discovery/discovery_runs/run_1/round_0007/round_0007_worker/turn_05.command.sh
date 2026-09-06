set -euo pipefail
echo 'shared lib tree:'
find /shared -maxdepth 3 -type f | sed -n '1,80p'
echo 'python sys.path:'
python - <<'PY'
import sys, os
print(sys.path)
print('PYTHONPATH', os.environ.get('PYTHONPATH'))
PY
