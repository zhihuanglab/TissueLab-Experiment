set -euo pipefail
python - <<'PY'
import sys, json
from pathlib import Path
sys.path.insert(0,'/shared/lib')
from shared_analysis import load_training_cohort, build_cell_table, load_region_polygons
cohort = load_training_cohort('/data')
row = cohort.iloc[0]
slide_path = Path('/data')/row['slide_name']
print('sample donor', row['donor_id'], row['slide_name'])
polys = load_region_polygons(slide_path, scale=16.0)
print('regions:', sorted(polys)[:20], '... total', len(polys))
print('CA1 polygon count:', len(polys.get('CA1', [])))
cells = build_cell_table(slide_path, include_regions=True, include_geometry=True)
print('columns:', list(cells.columns))
ca1 = cells[cells['region']=='CA1']
print('CA1 n=', len(ca1))
print('cell types top 12 in CA1:')
print(ca1['cell_type'].value_counts().head(12).to_string())
for ct in ['Reactive Astrocyte','Astrocyte']:
    s = ca1.loc[ca1['cell_type']==ct, 'area']
    print(ct, 'count', len(s), 'median area', float(s.median()) if len(s) else None)
PY
