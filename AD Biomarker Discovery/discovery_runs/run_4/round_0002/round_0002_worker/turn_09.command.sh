python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis.stats import bootstrap_partial_correlation
print(inspect.signature(bootstrap_partial_correlation))
print(inspect.getsource(bootstrap_partial_correlation)[:2000])
PY
