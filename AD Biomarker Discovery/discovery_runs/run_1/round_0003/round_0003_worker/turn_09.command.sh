set -e
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
from shared_analysis import stats, artifacts
print('stats funcs:', [n for n in dir(stats) if not n.startswith('_')])
print('\npartial_correlation signature:', inspect.signature(stats.partial_correlation))
print(inspect.getsource(stats.partial_correlation))
print('\nleave_one_out_summary signature:', inspect.signature(stats.leave_one_out_summary))
print(inspect.getsource(stats.leave_one_out_summary))
print('\nwrite_donor_feature_table signature:', inspect.signature(artifacts.write_donor_feature_table))
PY
