set -euo pipefail
echo '--- __init__.py ---'
sed -n '1,200p' /shared/lib/shared_analysis/__init__.py
echo '--- ls -lb ---'
ls -lb /shared/lib/shared_analysis | sed -n '1,80p'
python - <<'PY'
import os, sys
print(os.path.exists('/shared/lib/shared_analysis/__init__.py'))
print(os.listdir('/shared/lib/shared_analysis')[:20])
import importlib.machinery
spec = importlib.machinery.PathFinder.find_spec('shared_analysis', ['/shared/lib'])
print(spec)
PY
