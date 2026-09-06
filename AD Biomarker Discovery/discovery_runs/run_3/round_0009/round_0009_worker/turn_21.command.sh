set -euo pipefail
python - <<'PY'
import sys
sys.path.append('/shared/lib')
from shared_analysis import load_training_cohort
from pathlib import Path
cohort=load_training_cohort(Path('/data'))
print(cohort.columns.tolist())
PY
