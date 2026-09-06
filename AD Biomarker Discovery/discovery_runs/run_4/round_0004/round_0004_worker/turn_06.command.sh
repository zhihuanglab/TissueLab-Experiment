set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0,'/shared/lib')
from pathlib import Path
import inspect
from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort, load_region_polygons
print('build_cell_table sig:', inspect.signature(build_cell_table))
cohort = load_training_cohort('/data')
slide = cohort.iloc[0]['slide_name']
slide_path = Path('/data')/slide
print('slide:', slide)
cells = build_cell_table(slide_path, include_regions=False, include_geometry=True)
print('n_cells:', len(cells))
print('cols:', list(cells.columns))
print(cells.head(2).to_dict(orient='records'))
polys = load_region_polygons(slide_path)
print('polygon regions:', list(polys)[:10])
for region, plist in polys.items():
    print('region', region, 'n_polygons', len(plist))
    if plist:
        arr = plist[0]
        print('first polygon shape', getattr(arr, 'shape', None), 'first rows', arr[:3])
    break
PY
