set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis import sea_ad_lfb
for name in ['compute_contour_geometry','build_cell_table','load_region_polygons']:
    obj=getattr(sea_ad_lfb,name)
    print(f'--- {name} ---')
    src=inspect.getsource(obj)
    print('\n'.join(src.splitlines()[:80]))
PY
