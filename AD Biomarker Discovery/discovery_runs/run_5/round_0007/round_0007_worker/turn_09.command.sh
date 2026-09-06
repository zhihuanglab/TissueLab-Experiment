set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
import shared_analysis.stats as st
print('members:', [m for m in dir(st) if not m.startswith('_')])
for name in ['residualized_loo_predictive_correlation','leave_one_out_partial_correlation']:
    print(f'\n{name} sig:', inspect.signature(getattr(st,name)))
    src=inspect.getsource(getattr(st,name))
    print(src[:5000])
PY
