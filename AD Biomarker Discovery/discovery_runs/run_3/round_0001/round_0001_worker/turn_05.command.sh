set -euo pipefail
echo 'shared dirs:'
find /shared -maxdepth 3 -type d | sed -n '1,80p'
echo 'shared files matching *analysis*'
find /shared -maxdepth 4 \( -type f -o -type d \) | grep -E 'shared_analysis|template|runtime_quickstart' | sed -n '1,120p'
python - <<'PY'
import sys, os
print('sys.path head:', sys.path[:10])
PY
