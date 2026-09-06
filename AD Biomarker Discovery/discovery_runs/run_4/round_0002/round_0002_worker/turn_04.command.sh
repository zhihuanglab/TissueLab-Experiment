set -e
echo 'PYTHONPATH='$PYTHONPATH
ls -la /shared
find /shared -maxdepth 2 -type f | sed -n '1,40p'
python - <<'PY'
import sys, os
print(sys.path)
PY
