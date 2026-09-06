set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0, '/shared/lib')
from pathlib import Path
from shared_analysis import build_cell_table
print('import ok')

for q in Path('/data').rglob('round_0009_worker/result.py'):
    print('FOUND PRIOR RESULT', q)
    print(q.read_text()[:2800])
    break
else:
    print('No local prior result under /data')

slide = Path('/data/H20.33.024-A12-LFB.svs.zarr')
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print('columns:', cells.columns.tolist())
print('shape:', cells.shape)
print(cells[['cell_label','region_label']].head().to_string())
print('region counts top:', cells['region_label'].value_counts(dropna=False).head(10).to_dict())
print('labels top:', cells['cell_label'].value_counts().head(10).to_dict())
PY
