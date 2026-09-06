set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis import build_cell_table
import shared_analysis.sea_ad_lfb as m
print(inspect.getsource(m.build_cell_table))
PY
