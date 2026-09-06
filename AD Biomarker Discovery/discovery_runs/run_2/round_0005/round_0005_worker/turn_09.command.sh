set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
from shared_analysis.stats import residualized_loo_predictive_correlation
print('\n'.join(inspect.getsource(residualized_loo_predictive_correlation).splitlines()[:220]))
PY
