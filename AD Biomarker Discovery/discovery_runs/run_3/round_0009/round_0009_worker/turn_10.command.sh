set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis import stats, artifacts
for name in ['partial_correlation','leave_one_out_summary','residualized_loo_predictive_correlation']:
    fn=getattr(stats,name)
    print('\n',name, inspect.signature(fn))
    try:
        src=inspect.getsource(fn)
        print(src[:2500])
    except Exception as e:
        print('source err',e)
for name in ['write_donor_feature_table','build_results_payload','write_results_payload']:
    fn=getattr(artifacts,name)
    print('\n',name, inspect.signature(fn))
PY
