set -euo pipefail
python - <<'PY'
import sys, os
print('PYTHONPATH env', os.environ.get('PYTHONPATH'))
sys.path.insert(0,'/shared/lib')
import shared_analysis
print(shared_analysis)
PY
