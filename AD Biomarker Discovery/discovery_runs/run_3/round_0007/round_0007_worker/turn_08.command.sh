set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
from shared_analysis.sea_ad_lfb import assign_centroids_to_regions
print(inspect.getsource(assign_centroids_to_regions))
PY
