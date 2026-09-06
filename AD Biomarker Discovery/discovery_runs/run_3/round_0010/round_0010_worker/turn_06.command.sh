set -euo pipefail
ls -R /shared | head -200
python - <<'PY'
import os, sys
print('PYTHONPATH env', os.environ.get('PYTHONPATH'))
for p in ['/shared/lib', '/shared', '/shared/templates']:
    print(p, os.path.exists(p), os.listdir(p)[:20] if os.path.isdir(p) else None)
PY
