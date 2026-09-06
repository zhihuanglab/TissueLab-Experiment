set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.append('/shared/lib')
from shared_analysis import artifacts
for name in ['write_donor_feature_table','build_results_payload','write_results_payload','validate_results_payload','coerce_results_payload']:
    fn=getattr(artifacts,name)
    print(f'\n## {name}{inspect.signature(fn)}')
    print(inspect.getsource(fn)[:2500])
PY
