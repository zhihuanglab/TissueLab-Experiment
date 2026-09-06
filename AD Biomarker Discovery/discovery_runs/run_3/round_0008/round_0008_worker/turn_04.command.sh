set -euo pipefail
python - <<'PY'
import sys, os
print(sys.path[:5])
print('shared exists', os.path.exists('/shared/lib'))
import os
for root,dirs,files in os.walk('/shared/lib'):
    print(root, files[:5])
    if root.count('/')>3:
        break
PY
