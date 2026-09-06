set -e
python - <<'PY'
import sys
sys.path.insert(0, '/shared/lib')
from pathlib import Path
from shared_analysis.sea_ad_lfb import build_cell_table, load_region_polygons
slide = Path('/data/H19.33.004-A12-LFB.svs.zarr')
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print('cell columns:', list(cells.columns))
print(cells[['x_centroid','y_centroid','cell_type','region_name']].head(3).to_string())
print('region counts top10:')
print(cells['region_name'].value_counts(dropna=False).head(10).to_string())
polys = load_region_polygons(slide)
print('\npolygon region keys:', list(polys.keys())[:10])
ca1 = polys.get('CA1')
print('CA1 polygons n=', None if ca1 is None else len(ca1))
if ca1:
    print('first CA1 point=', ca1[0][:3])
PY
