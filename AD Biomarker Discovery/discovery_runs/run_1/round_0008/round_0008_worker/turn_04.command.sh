set -euo pipefail
python - <<'PY'
import sys, os
print(sys.path)
print(os.listdir('/shared'))
print(os.listdir('/shared/lib')[:20])
PY
