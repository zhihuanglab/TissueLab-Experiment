set -e
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis import stats, artifacts
for name in ['partial_correlation','residualized_loo_predictive_correlation']:
    fn=getattr(stats,name)
    print(f'\n{name}{inspect.signature(fn)}')
    print(inspect.getsource(fn)[:2500])
for name in ['write_donor_feature_table','build_results_payload','validate_results_payload']:
    fn=getattr(artifacts,name)
    print(f'\n{name}{inspect.signature(fn)}')
    print(inspect.getsource(fn)[:1800])
PY
