set -euo pipefail
python - <<'PY'
import sys, inspect, textwrap
sys.path.insert(0, '/shared/lib')
import shared_analysis.stats as st
for name in ['residualized_loo_predictive_correlation','residualize','leave_one_out_partial_correlation','bootstrap_partial_correlation']:
    print(f'--- {name} ---')
    print(textwrap.dedent(inspect.getsource(getattr(st,name)))[:5000])
PY
