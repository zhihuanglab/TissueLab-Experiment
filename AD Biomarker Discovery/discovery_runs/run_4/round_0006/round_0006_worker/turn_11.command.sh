set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis.sea_ad_lfb import load_region_annotations, open_slide_zarr
print(inspect.getsource(load_region_annotations)[:2000])
print(inspect.getsource(open_slide_zarr))
PY
