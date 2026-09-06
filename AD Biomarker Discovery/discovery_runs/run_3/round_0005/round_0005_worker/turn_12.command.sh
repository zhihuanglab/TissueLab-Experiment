set -e
python - <<'PY'
from pathlib import Path
import pandas as pd, zarr, json
cohort = pd.read_csv('/data/training_cohort.csv')
slide = Path('/data') / cohort.iloc[0]['slide_name']
root = zarr.open_group(str(slide), mode='r')
s = root['CustomAnnotations']['CA1_01']['annotation_json'][()]
print(type(s), len(str(s)))
print(str(s)[:1500])
PY
