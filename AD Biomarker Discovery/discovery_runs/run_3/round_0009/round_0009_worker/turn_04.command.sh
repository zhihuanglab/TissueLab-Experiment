set -euo pipefail
python - <<'PY'
import sys, os
print(sys.path[:10])
print('shared exists', os.path.exists('/shared/lib'))
import os
for root, dirs, files in os.walk('/shared/lib'):
    if root.count(os.sep) < 4:
        print(root, files[:5])
PY
