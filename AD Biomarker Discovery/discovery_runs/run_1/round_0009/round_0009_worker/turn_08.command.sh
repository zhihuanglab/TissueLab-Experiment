set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
import shared_analysis
print(inspect.getsource(shared_analysis))
from shared_analysis import stats, artifacts
print('--- stats functions ---')
for name in ['partial_correlation','leave_one_out_summary','residualized_loo_predictive_correlation']:
    fn=getattr(stats,name)
    print(f'\n{name}{inspect.signature(fn)}')
    print(inspect.getsource(fn)[:3000])
print('\n--- artifacts functions ---')
for name in ['write_donor_feature_table','build_results_payload','write_results_payload','validate_results_payload']:
    fn=getattr(artifacts,name)
    print(f'\n{name}{inspect.signature(fn)}')
    print(inspect.getsource(fn)[:2000])
PY
