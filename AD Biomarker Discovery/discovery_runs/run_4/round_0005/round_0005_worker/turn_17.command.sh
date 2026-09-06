set -euo pipefail
python - <<'PY'
import pandas as pd
cohort=pd.read_csv('/data/training_cohort.csv')
print(cohort.columns.tolist())
PY
