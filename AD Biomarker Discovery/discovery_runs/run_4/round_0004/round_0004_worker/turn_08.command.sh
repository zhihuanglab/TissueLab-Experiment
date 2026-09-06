set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis import stats, artifacts
for name in ['partial_correlation','leave_one_out_summary','residualized_loo_predictive_correlation']:
    print(f'--- {name} ---')
    print('\n'.join(inspect.getsource(getattr(stats,name)).splitlines()[:120]))
for name in ['write_donor_feature_table','build_results_payload','write_results_payload','validate_results_payload']:
    print(f'--- {name} ---')
    print('\n'.join(inspect.getsource(getattr(artifacts,name)).splitlines()[:120]))
PY
