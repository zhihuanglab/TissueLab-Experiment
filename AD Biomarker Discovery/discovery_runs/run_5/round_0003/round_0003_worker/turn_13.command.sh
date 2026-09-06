python - <<'PY'
import pandas as pd
cohort=pd.read_csv('/data/training_cohort.csv')
print(cohort.columns.tolist())
print(cohort.head(2).to_string())
PY
