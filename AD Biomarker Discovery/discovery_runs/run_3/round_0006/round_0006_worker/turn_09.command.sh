set -euo pipefail
python - <<'PY'
import sys, inspect, textwrap
sys.path.insert(0,'/shared/lib')
from shared_analysis import artifacts
for name in ['write_donor_feature_table','write_results_payload','build_results_payload','validate_results_payload']:
    fn=getattr(artifacts,name)
    print('---', name, inspect.signature(fn))
    print(textwrap.dedent('\n'.join(inspect.getsource(fn).splitlines()[:120])))
    print()
PY
