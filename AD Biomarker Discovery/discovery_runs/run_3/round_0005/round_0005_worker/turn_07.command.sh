set -e
python - <<'PY'
from pathlib import Path
import pandas as pd, zarr
cohort = pd.read_csv('/data/training_cohort.csv')
slide = Path('/data') / cohort.iloc[0]['slide_name']
root = zarr.open_group(str(slide), mode='r')
print(root.tree())
for name in ['SegmentationNode','ClassificationNode','CustomAnnotations']:
    if name in root:
        grp = root[name]
        print("\n==", name, "==")
        try:
            print(grp.tree())
        except Exception as e:
            print("tree failed", e)
PY
