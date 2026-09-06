set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis import leave_one_out_summary, write_donor_feature_table
print(inspect.signature(leave_one_out_summary))
print(inspect.getsource(leave_one_out_summary)[:5000])
print(inspect.signature(write_donor_feature_table))
PY
