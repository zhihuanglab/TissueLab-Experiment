set -euo pipefail
python - <<'PY'
import sys, inspect, textwrap
sys.path.insert(0,'/shared/lib')
from shared_analysis.sea_ad_lfb import load_training_cohort
print(textwrap.dedent('\n'.join(inspect.getsource(load_training_cohort).splitlines()[:120])))
PY
