set -e
python - <<'PY'
import pandas as pd, json
from pathlib import Path
cohort = pd.read_csv('/data/training_cohort.csv')
print('columns:', list(cohort.columns))
print('results exists:', Path('/scratch/results.json').exists())
print('report exists:', Path('/scratch/report.md').exists())
print('run_stdout exists:', Path('/scratch/run_stdout.txt').exists())
PY
