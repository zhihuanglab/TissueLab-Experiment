set -euo pipefail
python - <<'PY'
import sys, inspect, textwrap
sys.path.insert(0,'/shared/lib')
from shared_analysis.sea_ad_lfb import load_class_lookup
print(textwrap.dedent(inspect.getsource(load_class_lookup)))
PY
