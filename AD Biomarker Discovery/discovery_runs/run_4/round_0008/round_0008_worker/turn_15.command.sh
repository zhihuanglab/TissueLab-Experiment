set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
import shared_analysis.sea_ad_lfb as sea
for name in ['load_centroids','load_class_ids','load_class_lookup']:
    print(name, inspect.signature(getattr(sea,name)))
PY
