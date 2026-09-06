set -euo pipefail
python - <<'PY'
import inspect
from shared_analysis import build_cell_table, compute_contour_geometry
from shared_analysis import partial_correlation, residualized_loo_predictive_correlation, write_donor_feature_table
import shared_analysis
import shared_analysis.sea_ad_lfb as lfb
import shared_analysis.stats as stats
import shared_analysis.artifacts as artifacts

for name, obj in [('build_cell_table', lfb.build_cell_table),
                  ('compute_contour_geometry', lfb.compute_contour_geometry),
                  ('partial_correlation', stats.partial_correlation),
                  ('residualized_loo_predictive_correlation', stats.residualized_loo_predictive_correlation),
                  ('write_donor_feature_table', artifacts.write_donor_feature_table)]:
    print(f'\n=== {name} signature ===')
    print(inspect.signature(obj))
    try:
        src = inspect.getsource(obj)
        print('\n'.join(src.splitlines()[:120]))
    except Exception as e:
        print('source unavailable', e)
PY
