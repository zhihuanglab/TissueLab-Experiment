set -e
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis.stats import partial_correlation, bootstrap_partial_correlation_stability
print(inspect.signature(partial_correlation))
print(inspect.getsource(partial_correlation)[:1500])
print(inspect.signature(bootstrap_partial_correlation_stability))
print(inspect.getsource(bootstrap_partial_correlation_stability)[:2000])
PY
