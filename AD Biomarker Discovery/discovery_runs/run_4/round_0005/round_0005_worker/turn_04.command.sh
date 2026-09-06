set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
python - <<'PY'
import inspect, textwrap
from shared_analysis.sea_ad_lfb import build_cell_table, load_region_polygons, load_training_cohort
from shared_analysis.stats import partial_correlation, leave_one_out_summary, residualized_loo_predictive_correlation
from shared_analysis.artifacts import write_donor_feature_table, write_results_payload, build_results_payload, validate_results_payload
for fn in [load_training_cohort, build_cell_table, load_region_polygons, partial_correlation, leave_one_out_summary, residualized_loo_predictive_correlation, write_donor_feature_table, build_results_payload, validate_results_payload]:
    print(f"\n### SOURCE: {fn.__module__}.{fn.__name__}")
    src = inspect.getsource(fn)
    print(textwrap.dedent(src[:3500]))
PY
