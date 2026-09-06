set -euo pipefail
python - <<'PY'
import sys, inspect, textwrap
sys.path.insert(0,'/shared/lib')
from shared_analysis import sea_ad_lfb, stats
print('--- compute_contour_geometry ---')
print(textwrap.dedent(inspect.getsource(sea_ad_lfb.compute_contour_geometry)))
print('--- build_cell_table excerpt ---')
src=inspect.getsource(sea_ad_lfb.build_cell_table)
print(textwrap.dedent('\n'.join(src.splitlines()[:120])))
print('--- stats partial/loo signatures ---')
for name in ['partial_correlation','leave_one_out_summary','residualized_loo_predictive_correlation']:
    fn=getattr(stats,name)
    print(name, inspect.signature(fn))
    print(textwrap.dedent('\n'.join(inspect.getsource(fn).splitlines()[:80])))
    print()
PY
