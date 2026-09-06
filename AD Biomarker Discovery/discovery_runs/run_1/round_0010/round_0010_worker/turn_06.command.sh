set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.append('/shared/lib')
from shared_analysis import sea_ad_lfb, stats
for name in ['load_training_cohort','build_cell_table','compute_contour_geometry','load_region_polygons']:
    fn=getattr(sea_ad_lfb,name)
    print(f'\n## {name}{inspect.signature(fn)}')
    print(inspect.getsource(fn)[:2500])
for name in ['partial_correlation','leave_one_out_summary','residualized_loo_predictive_correlation']:
    fn=getattr(stats,name)
    print(f'\n## {name}{inspect.signature(fn)}')
    print(inspect.getsource(fn)[:2200])
PY
