set -euo pipefail
python - <<'PY'
import sys, os, glob
print(sys.path)
print('shared dirs', glob.glob('/shared/lib/*')[:20])
PY
