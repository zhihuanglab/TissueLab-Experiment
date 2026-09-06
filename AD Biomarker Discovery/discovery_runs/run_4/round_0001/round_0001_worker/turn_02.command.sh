set -euo pipefail
python - <<'PY'
import json
with open('/scratch/context_bundle.json') as f:
    obj = json.load(f)
print("TOP-LEVEL KEYS:", list(obj))
for k,v in obj.items():
    if isinstance(v, (dict,list)):
        s=json.dumps(v, indent=2)
        print("\n##", k)
        print(s[:20000])
    else:
        print("\n##", k, "=", v)
PY
python - <<'PY'
from pathlib import Path
import pandas as pd
p=Path('/data/training_cohort.csv')
print("\ntraining_cohort exists:", p.exists())
if p.exists():
    df=pd.read_csv(p)
    print(df.head().to_string())
    print(df.columns.tolist())
    print(df.shape)
PY
