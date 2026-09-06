set -euo pipefail
echo '--- /shared ---'
find /shared -maxdepth 3 -type f | sed -n '1,80p'
echo '--- python path probe ---'
python - <<'PY'
import sys, os
print(sys.path)
print('PYTHONPATH', os.environ.get('PYTHONPATH'))
PY
