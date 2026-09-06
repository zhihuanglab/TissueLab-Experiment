set -euo pipefail
python - <<'PY'
import sys, inspect, textwrap
sys.path.insert(0, '/shared/lib')
import shared_analysis.sea_ad_lfb as sea
print(textwrap.dedent(inspect.getsource(sea.assign_centroids_to_regions))[:3000])
PY
