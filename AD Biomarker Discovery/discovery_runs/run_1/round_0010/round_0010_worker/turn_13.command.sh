set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.append('/shared/lib')
from shared_analysis import stats
for name in ['residualize','_design_matrix','_clean_frame']:
    fn=getattr(stats,name)
    print(f'\n## {name}{inspect.signature(fn)}')
    print(inspect.getsource(fn)[:2500])
PY
