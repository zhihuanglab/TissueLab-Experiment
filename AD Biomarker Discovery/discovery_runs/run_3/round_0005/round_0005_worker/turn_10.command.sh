set -e
python - <<'PY'
from pathlib import Path
import pandas as pd, zarr
cohort = pd.read_csv('/data/training_cohort.csv')
slide = Path('/data') / cohort.iloc[0]['slide_name']
root = zarr.open_group(str(slide), mode='r')
for ann in ['CA1_01','CA2_01']:
    grp=root['CustomAnnotations'][ann]
    print(ann, list(grp.keys()))
    for k in grp.keys():
        a=grp[k]
        print(" ",k, getattr(a,'shape',None), getattr(a,'dtype',None))
        try:
            print("   head", a[:5].tolist())
        except Exception:
            pass
PY
