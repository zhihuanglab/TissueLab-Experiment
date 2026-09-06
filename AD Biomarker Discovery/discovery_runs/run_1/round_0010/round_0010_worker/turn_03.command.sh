set -euo pipefail
echo '--- shared_analysis.sea_ad_lfb signatures/head ---'
python - <<'PY'
import inspect
from shared_analysis import sea_ad_lfb
for name in ['load_training_cohort','build_cell_table','compute_contour_geometry','load_region_polygons']:
    fn=getattr(sea_ad_lfb,name)
    print(f'\n## {name}{inspect.signature(fn)}')
    src=inspect.getsource(fn)
    print(src[:4000])
PY
echo
echo '--- stats helpers signatures/head ---'
python - <<'PY'
import inspect
from shared_analysis import stats
for name in ['partial_correlation','leave_one_out_summary','residualized_loo_predictive_correlation']:
    fn=getattr(stats,name)
    print(f'\n## {name}{inspect.signature(fn)}')
    src=inspect.getsource(fn)
    print(src[:3000])
PY
