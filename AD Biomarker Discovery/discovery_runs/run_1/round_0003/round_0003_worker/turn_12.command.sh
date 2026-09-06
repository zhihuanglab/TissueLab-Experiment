set -e
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
from shared_analysis.artifacts import write_donor_feature_table
print(inspect.getsource(write_donor_feature_table))
PY
