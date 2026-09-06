set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
from shared_analysis.sea_ad_lfb import load_region_polygons
print(inspect.getsource(load_region_polygons))
PY
