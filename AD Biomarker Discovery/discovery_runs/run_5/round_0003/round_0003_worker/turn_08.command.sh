set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.append('/shared/lib')
from shared_analysis import stats
for name in ['partial_correlation','residualized_loo_predictive_correlation','leave_one_out_summary']:
    fn=getattr(stats,name)
    print('---', name, inspect.signature(fn), '---')
    src=inspect.getsource(fn).splitlines()
    for line in src[:220]:
        print(line)
    print()
PY
