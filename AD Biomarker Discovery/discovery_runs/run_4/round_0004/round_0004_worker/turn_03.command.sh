set -euo pipefail
python - <<'PY'
from pathlib import Path
import json
from shared_analysis import build_cell_table, load_training_cohort, load_region_polygons
cohort = load_training_cohort('/data')
slide = cohort.iloc[0]['slide_name']
slide_path = Path('/data')/slide
print('slide:', slide)
cells = build_cell_table(slide_path, include_regions=False, include_geometry=True)
print('cell columns:', list(cells.columns)[:40])
print('n_cells:', len(cells))
print(cells.head(2).to_string())
polys = load_region_polygons(slide_path)
print('polygon keys sample:', list(polys.keys())[:10])
for k,v in polys.items():
    print('first polygon key', k, 'type', type(v), 'len', len(v) if hasattr(v,'__len__') else 'na')
    print('sample first item type', type(v[0]) if len(v) else None)
    break
PY
