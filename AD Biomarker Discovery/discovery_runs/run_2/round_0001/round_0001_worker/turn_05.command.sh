python - <<'PY'
import sys, os
print('PYTHONPATH', os.environ.get('PYTHONPATH'))
print('sys.path head', sys.path[:8])
import pathlib
print('/shared/lib exists', pathlib.Path('/shared/lib').exists())
PY

echo 'shared/lib entries:'
find /shared/lib -maxdepth 2 -type f | sed -n '1,40p'
