set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis import sea_ad_lfb
print('build_cell_table source excerpt:')
src = inspect.getsource(sea_ad_lfb.build_cell_table)
print(src[:4000])
print('\ncompute_contour_geometry source excerpt:')
print(inspect.getsource(sea_ad_lfb.compute_contour_geometry)[:4000])
PY
