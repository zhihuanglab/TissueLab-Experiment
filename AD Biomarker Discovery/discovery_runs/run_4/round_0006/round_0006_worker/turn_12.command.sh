set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis import partial_correlation, residualized_loo_predictive_correlation, leave_one_out_summary
print(inspect.signature(partial_correlation))
print(inspect.getsource(partial_correlation)[:3000])
print(inspect.signature(residualized_loo_predictive_correlation))
print(inspect.getsource(residualized_loo_predictive_correlation)[:4000])
PY
