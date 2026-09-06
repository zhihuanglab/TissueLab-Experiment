set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
import shared_analysis.sea_ad_lfb as lfb
import shared_analysis.stats as stats
import shared_analysis.artifacts as artifacts

targets = [
    ('build_cell_table', lfb.build_cell_table),
    ('compute_contour_geometry', lfb.compute_contour_geometry),
    ('partial_correlation', stats.partial_correlation),
    ('residualized_loo_predictive_correlation', stats.residualized_loo_predictive_correlation),
    ('leave_one_out_summary', stats.leave_one_out_summary),
    ('write_donor_feature_table', artifacts.write_donor_feature_table),
]
for name, obj in targets:
    print(f'\n=== {name} signature ===')
    print(inspect.signature(obj))
    src = inspect.getsource(obj).splitlines()
    for line in src[:120]:
        print(line)
PY
