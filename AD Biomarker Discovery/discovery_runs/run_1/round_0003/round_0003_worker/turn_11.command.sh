set -e
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
from shared_analysis.sea_ad_lfb import build_cell_table
print(inspect.getsource(build_cell_table)[:4000])
PY
