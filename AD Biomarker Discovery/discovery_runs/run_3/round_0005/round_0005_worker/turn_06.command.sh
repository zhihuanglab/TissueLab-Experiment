set -e
python - <<'PY'
from pathlib import Path
import pandas as pd, json, os, zarr
cohort = pd.read_csv('/data/training_cohort.csv')
slide = Path('/data') / cohort.iloc[0]['slide_name']
print("slide:", slide)
root = zarr.open(str(slide), mode='r')
def walk(group, prefix=''):
    for k,v in group.items():
        path=f"{prefix}/{k}" if prefix else k
        if hasattr(v,'shape'):
            print(path, 'shape=', v.shape, 'dtype=', v.dtype)
        else:
            print(path, '(group)')
            walk(v, path)
walk(root)
PY
