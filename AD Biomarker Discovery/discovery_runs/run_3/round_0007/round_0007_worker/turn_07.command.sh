set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
from shared_analysis.stats import residualized_loo_predictive_correlation
print(inspect.getsource(residualized_loo_predictive_correlation))
PY
