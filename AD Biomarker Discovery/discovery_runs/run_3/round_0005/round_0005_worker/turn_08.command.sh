set -e
python - <<'PY'
from pathlib import Path
import pandas as pd, zarr
cohort = pd.read_csv('/data/training_cohort.csv')
slide = Path('/data') / cohort.iloc[0]['slide_name']
root = zarr.open_group(str(slide), mode='r')
print(type(root))
print([m for m in dir(root) if not m.startswith('_')][:80])
print("keys", list(root.keys()))
for k in root.keys():
    print(k, type(root[k]))
    grp=root[k]
    try:
        print(" subkeys", list(grp.keys())[:20])
    except Exception as e:
        print(" no subkeys", e)
PY
