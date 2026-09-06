set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
python - <<'PY'
from pathlib import Path
from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort, load_region_polygons
cohort = load_training_cohort('/data')
slide = cohort.iloc[0]['slide_name']
slide_path = Path('/data')/slide
print('slide:', slide)
cells = build_cell_table(slide_path, include_regions=False, include_geometry=True)
print('n_cells:', len(cells))
print('cell columns:', list(cells.columns))
print(cells.head(2).to_string())
polys = load_region_polygons(slide_path)
print('polygon keys sample:', list(polys.keys())[:10])
for k,v in polys.items():
    print('first polygon key', k, 'type', type(v), 'len', len(v) if hasattr(v,'__len__') else 'na')
    print('sample first item type', type(v[0]) if len(v) else None)
    if len(v):
        print('sample first polygon first 3 pts:', v[0][:3])
    break
PY
