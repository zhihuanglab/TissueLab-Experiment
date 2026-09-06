set -euo pipefail
echo "PYTHONPATH=$PYTHONPATH"
ls -la /shared
find /shared -maxdepth 3 -type d | sed -n '1,40p'
python - <<'PY'
import sys
print(sys.path[:10])
PY
