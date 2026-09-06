set -e
python - <<'PY'
import sys
sys.path.insert(0,'/shared/lib')
from shared_analysis import load_training_cohort
c=load_training_cohort('/data')
print(list(c.columns))
PY
