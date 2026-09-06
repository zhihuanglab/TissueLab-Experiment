set -euo pipefail
python - <<'PY'
import sys, json
from pathlib import Path
sys.path.insert(0, '/shared/lib')
from shared_analysis import build_cell_table

slide = Path('/data/H20.33.024-A12-LFB.svs.zarr')
cells = build_cell_table(slide, include_regions=True, include_geometry=True)
print('columns:', cells.columns.tolist())
print('n_cells:', len(cells))
print('regions counts top:', cells['region'].value_counts(dropna=False).head(10).to_dict())
print('cell types:', cells['cell_type'].value_counts().to_dict())
ca1 = cells[cells['region']=='CA1']
print('CA1 n:', len(ca1))
print('CA1 types:', ca1['cell_type'].value_counts().to_dict())
pyr = ca1[ca1['cell_type']=='Pyramidal Neuron']
ra = ca1[ca1['cell_type']=='Reactive Astrocyte']
print('CA1 pyramidal n:', len(pyr), 'area quantiles:', pyr['area'].quantile([0.25,0.33,0.5]).to_dict())
print('CA1 reactive astrocyte n:', len(ra))
print('head:', ca1[['x','y','cell_type','region','area']].head().to_dict(orient='records'))
PY
