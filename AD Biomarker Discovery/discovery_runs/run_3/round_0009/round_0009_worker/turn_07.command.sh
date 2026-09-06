set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis import sea_ad_lfb
print('load_training_cohort sig', inspect.signature(sea_ad_lfb.load_training_cohort))
print('build_cell_table sig', inspect.signature(sea_ad_lfb.build_cell_table))
src = inspect.getsource(sea_ad_lfb.build_cell_table)
print(src[:4000])
PY
