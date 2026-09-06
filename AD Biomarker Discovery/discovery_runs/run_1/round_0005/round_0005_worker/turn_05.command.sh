set -e
echo 'shared listing'
find /shared -maxdepth 3 -type f | sed -n '1,120p'
echo 'python paths check'
python - <<'PY'
import sys, os
print(sys.path[:10])
print('shared exists', os.path.exists('/shared'))
print('shared/lib exists', os.path.exists('/shared/lib'))
PY
