set -e
python - <<'PY'
from shared_analysis import build_cell_table
from pathlib import Path
import inspect
print('build_cell_table signature:', inspect.signature(build_cell_table))
src = inspect.getsource(build_cell_table)
for needle in ['include_regions', 'region', '16', 'CustomAnnotations']:
    print(f'contains {needle!r}:', needle in src)
print(src[:5000])
PY
echo
python - <<'PY'
from shared_analysis import build_cell_table
from pathlib import Path
slide = Path('/data/H19.33.004-A12-LFB.svs.zarr')
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print('columns:', list(cells.columns))
print('n=', len(cells), 'region nonnull=', int(cells['region'].notna().sum()))
print(cells[['x','y','cell_type','region']].head().to_string())
print('cell types:', cells['cell_type'].value_counts().head(10).to_dict())
print('regions:', cells['region'].value_counts(dropna=False).to_dict())
PY
