set -euo pipefail
python - <<'PY'
import sys, os
print(sys.path)
print('shared exists', os.path.exists('/shared'))
for p in ['/shared/lib','/shared/templates']:
    print(p, os.path.exists(p), os.listdir(p)[:10] if os.path.exists(p) else None)
PY
