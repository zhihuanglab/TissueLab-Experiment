set -euo pipefail
python - <<'PY'
import sys, inspect, re
sys.path.append('/shared/lib')
from shared_analysis import stats
import shared_analysis
import pkgutil, pathlib, os, subprocess, textwrap
print('residualized_loo_predictive_correlation sig:', inspect.signature(stats.residualized_loo_predictive_correlation))
print(inspect.getsource(stats.residualized_loo_predictive_correlation))
PY
grep -R "Adjusted Score" -n /shared/lib /shared/templates /scratch || true
grep -R "selection_score" -n /shared/lib /shared/templates | head -50 || true
