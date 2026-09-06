set -euo pipefail
python - <<'PY'
from pathlib import Path
import pandas as pd
from shared_analysis import build_cell_table

# inspect prior round artifacts if locally available
for p in [
    Path('/data/autoresearch_runs'),
    Path('/data')
]:
    print('exists', p, p.exists())

# find prior round 9 result if present
for q in Path('/data').rglob('round_0009_worker/result.py'):
    print('FOUND PRIOR RESULT', q)
    print(q.read_text()[:3500])
    break
else:
    print('No local prior round_0009_worker/result.py found under /data')

slide = Path('/data/H20.33.024-A12-LFB.svs.zarr')
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print('columns:', cells.columns.tolist())
print('shape:', cells.shape)
print(cells[['cell_label','region_label']].head())
print('region counts top:', cells['region_label'].value_counts(dropna=False).head(10).to_dict())
print('labels:', cells['cell_label'].value_counts().to_dict())
PY
