set -euo pipefail
PYTHONPATH=/shared/lib python - <<'PY'
import sys, os
print('PYTHONPATH env', os.environ.get('PYTHONPATH'))
print(sys.path[:5])
import shared_analysis
print(shared_analysis)
from shared_analysis import sea_ad_lfb
import inspect
print(inspect.signature(sea_ad_lfb.build_cell_table))
PY
