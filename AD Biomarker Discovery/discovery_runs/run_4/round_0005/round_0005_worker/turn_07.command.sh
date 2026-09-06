set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0,'/shared/lib')
from pathlib import Path
from shared_analysis.sea_ad_lfb import build_cell_table, load_region_polygons
slide=Path('/data/H20.33.024-A12-LFB.svs.zarr')
cells=build_cell_table(slide, include_regions=True, include_geometry=False)
print(cells.head())
print(cells.columns.tolist())
print(cells['region'].value_counts(dropna=False).head(12))
print(cells['cell_type'].value_counts().to_dict())
polys=load_region_polygons(slide)
print('regions', polys.keys())
print('CA1 polygons', len(polys.get('CA1', [])), 'first shape', polys.get('CA1',[None])[0].shape if polys.get('CA1') else None)
PY
