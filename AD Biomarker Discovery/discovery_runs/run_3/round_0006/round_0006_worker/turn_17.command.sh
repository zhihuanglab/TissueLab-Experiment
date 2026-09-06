set -euo pipefail
python - <<'PY'
import sys, inspect, textwrap
sys.path.insert(0,'/shared/lib')
from shared_analysis import sea_ad_lfb
for name in ['load_contours','load_centroids','load_class_ids','load_region_polygons','assign_centroids_to_regions','load_region_annotations']:
    fn=getattr(sea_ad_lfb,name)
    print('---', name, inspect.signature(fn))
    print(textwrap.dedent('\n'.join(inspect.getsource(fn).splitlines()[:120])))
    print()
PY
