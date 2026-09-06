set -e
find /shared -maxdepth 2 -type f ! -name '._*' | sed -n '1,80p'
python - <<'PY'
import sys, os
print('sys.path first entries:', sys.path[:5])
sys.path.insert(0,'/shared/lib')
from shared_analysis import sea_ad_lfb
import inspect
print('loaded shared_analysis from', sea_ad_lfb.__file__)
for name in ['build_cell_table','load_region_polygons']:
    fn = getattr(sea_ad_lfb, name)
    print(f'\n{name}{inspect.signature(fn)}')
    src = inspect.getsource(fn)
    print(src[:2200])
PY
