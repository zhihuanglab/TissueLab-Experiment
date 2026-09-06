set -e
python - <<'PY'
from pathlib import Path
import inspect
from shared_analysis.sea_ad_lfb import build_cell_table, load_region_polygons
slide = Path('/data/H19.33.004-A12-LFB.svs.zarr')
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print('cell columns:', list(cells.columns))
print('head:')
print(cells.head(3).to_string())
print('region counts top10:')
print(cells['region_name'].value_counts(dropna=False).head(10).to_string())
print('\nload_region_polygons sample keys:')
polys = load_region_polygons(slide)
print(type(polys), list(polys.keys())[:10])
ca1 = polys.get('CA1')
if ca1 is not None:
    print('CA1 polygons n=', len(ca1))
    first = ca1[0]
    print('first CA1 polygon shape=', getattr(first, 'shape', None), 'first points=', first[:3])
PY
