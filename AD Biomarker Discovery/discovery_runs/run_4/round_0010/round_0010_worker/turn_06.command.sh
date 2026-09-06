set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis import stats, artifacts
for mod, names in [(stats,['partial_correlation','leave_one_out_summary','residualized_loo_predictive_correlation']),
                   (artifacts,['build_results_payload','write_donor_feature_table','write_results_payload','coerce_results_payload','validate_results_payload'])]:
    print(f'\nMODULE {mod.__name__}')
    for name in names:
        fn=getattr(mod,name)
        print(f'\n{name}{inspect.signature(fn)}')
        src=inspect.getsource(fn)
        print('\n'.join(src.splitlines()[:60]))
PY
