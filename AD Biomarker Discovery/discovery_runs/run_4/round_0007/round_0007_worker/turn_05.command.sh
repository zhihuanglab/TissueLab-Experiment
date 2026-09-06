set -euo pipefail
echo 'shared lib contents'
find /shared/lib -maxdepth 2 -type f | sed -n '1,120p'
echo 'python path test'
python - <<'PY'
import sys, os
print(sys.version)
print('PYTHONPATH env=', os.environ.get('PYTHONPATH'))
print(sys.path[:10])
PY
echo 'scratch files'
find /scratch -maxdepth 2 -type f | sed -n '1,120p'
