set -euo pipefail
python - <<'PY'
import sys, inspect, textwrap
sys.path.insert(0,'/shared/lib')
from shared_analysis.stats import bootstrap_partial_correlation
print(inspect.signature(bootstrap_partial_correlation))
print(textwrap.dedent('\n'.join(inspect.getsource(bootstrap_partial_correlation).splitlines()[:140])))
PY
