set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
from shared_analysis.stats import partial_correlation, leave_one_out_summary, residualized_loo_predictive_correlation
from shared_analysis.artifacts import write_donor_feature_table
print('partial_correlation:\n', '\n'.join(inspect.getsource(partial_correlation).splitlines()[:80]))
print('\nleave_one_out_summary:\n', '\n'.join(inspect.getsource(leave_one_out_summary).splitlines()[:140]))
print('\nwrite_donor_feature_table:\n', '\n'.join(inspect.getsource(write_donor_feature_table).splitlines()[:120]))
PY
