set -euo pipefail
python - <<'PY'
from pathlib import Path
import inspect
from shared_analysis import build_cell_table, load_training_cohort
from shared_analysis.sea_ad_lfb import load_region_polygons
import pandas as pd

data_root = Path('/data')
cohort = load_training_cohort(data_root)
row = cohort.iloc[0]
slide_path = data_root / row['slide_name']
print('sample slide:', row['slide_name'])

cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
print('cell table shape:', cells.shape)
print('columns:', list(cells.columns))
print('class counts head:')
print(cells['classification'].value_counts().head(10).to_string())
if 'region_label' in cells.columns:
    print('region labels:', sorted(map(str, pd.Series(cells['region_label']).dropna().unique())))

print('\nload_region_polygons signature:')
print(inspect.signature(load_region_polygons))
polys = load_region_polygons(slide_path)
print('polygon keys:', list(polys.keys())[:10])
for k,v in list(polys.items())[:3]:
    print('first polygon example', k, 'n=', len(v), 'first shape len=', len(v[0]) if v else None)
PY
