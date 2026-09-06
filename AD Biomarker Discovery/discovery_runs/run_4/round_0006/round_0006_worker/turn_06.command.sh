set -euo pipefail
find /shared/lib -maxdepth 3 \( -type f -o -type d \) | sed -n '1,120p'
python - <<'PY'
import sys, importlib.util
for p in sys.path:
    print(p)
spec=importlib.util.find_spec('shared_analysis')
print('spec=',spec)
PY
