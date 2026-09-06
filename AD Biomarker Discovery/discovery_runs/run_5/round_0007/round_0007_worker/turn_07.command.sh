set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis import stats, artifacts, load_training_cohort
print('load_training_cohort source:')
print(inspect.getsource(load_training_cohort))
print('\npartial_correlation signature:', inspect.signature(stats.partial_correlation))
print(inspect.getsource(stats.partial_correlation))
print('\nleave_one_out_summary signature:', inspect.signature(stats.leave_one_out_summary))
print(inspect.getsource(stats.leave_one_out_summary))
print('\nwrite_donor_feature_table signature:', inspect.signature(artifacts.write_donor_feature_table))
print(inspect.getsource(artifacts.write_donor_feature_table))
PY
