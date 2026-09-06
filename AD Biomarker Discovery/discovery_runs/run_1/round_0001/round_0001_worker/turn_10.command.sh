set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0,'/shared/lib')
from shared_analysis import load_training_cohort, load_class_names
cohort=load_training_cohort('/data')
slide='/data/' + cohort.iloc[0]['slide_name']
print(load_class_names(slide))
PY
