set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
python - <<'PY'
import inspect
from shared_analysis import sea_ad_lfb, stats
for name in ['load_training_cohort','build_cell_table','compute_contour_geometry','load_region_polygons']:
    fn=getattr(sea_ad_lfb,name)
    print(f'\n## {name}{inspect.signature(fn)}')
    src=inspect.getsource(fn)
    print(src[:3000])
for name in ['partial_correlation','leave_one_out_summary','residualized_loo_predictive_correlation']:
    fn=getattr(stats,name)
    print(f'\n## {name}{inspect.signature(fn)}')
    src=inspect.getsource(fn)
    print(src[:2500])
PY
