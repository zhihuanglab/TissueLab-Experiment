set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.append('/shared/lib')
from shared_analysis import stats, artifacts
from shared_analysis import build_cell_table
print('partial_correlation signature:', inspect.signature(stats.partial_correlation))
print(inspect.getsource(stats.partial_correlation))
print('leave_one_out_summary signature:', inspect.signature(stats.leave_one_out_summary))
print(inspect.getsource(stats.leave_one_out_summary)[:4000])
print('build_results_payload signature:', inspect.signature(artifacts.build_results_payload))
print(inspect.getsource(artifacts.write_donor_feature_table)[:2500])
PY
