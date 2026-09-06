set -euo pipefail
python - <<'PY'
import sys, os
print(sys.path[:10])
PY
echo 'shared lib dirs:'
find /shared/lib -maxdepth 2 -type f | sed -n '1,80p'
echo 'template helper imports maybe package layout:'
python - <<'PY'
import os
for root, dirs, files in os.walk('/shared/lib'):
    if root.count(os.sep) > 4:
        continue
    for f in files:
        if f.endswith('.py') and ('shared_analysis' in root or f=='__init__.py'):
            print(os.path.join(root,f))
PY
