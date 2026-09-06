set -e
python - <<'PY'
import sys
sys.path.insert(0, '/shared/lib')
from shared_analysis import build_cell_table
from pathlib import Path
import inspect
print('build_cell_table signature:', inspect.signature(build_cell_table))
src = inspect.getsource(build_cell_table)
print('contains divide by 16 hint:', '/ 16' in src or '16.0' in src)
print(src[:3000])
print('--- one slide sample ---')
slide = Path('/data/H19.33.004-A12-LFB.svs.zarr')
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print('columns:', list(cells.columns))
print('n=', len(cells), 'region nonnull=', int(cells["region"].notna().sum()))
print(cells[['x','y','cell_type','region']].head().to_string())
print('regions sample:', cells['region'].value_counts(dropna=False).head(10).to_dict())
PY
