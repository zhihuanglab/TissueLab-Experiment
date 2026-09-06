set -euo pipefail
echo 'shared lib files:'
find /shared/lib -maxdepth 2 -type f | sed -n '1,120p'
echo 'python path and available dirs'
python - <<'PY'
import sys, os
print(sys.path[:10])
print(os.listdir('/shared/lib')[:20])
PY
