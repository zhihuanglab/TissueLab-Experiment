set -euo pipefail
python - <<'PY'
import sys, inspect, zarr
sys.path.insert(0, '/shared/lib')
from pathlib import Path
from shared_analysis.sea_ad_lfb import build_cell_table, load_region_polygons

print('--- build_cell_table source snippet ---')
src = inspect.getsource(build_cell_table)
print(src[:3500])

slide = Path('/data/H19.33.004-A12-LFB.svs.zarr')
print('\n--- sample cell table ---')
df = build_cell_table(slide, include_regions=True, include_geometry=False)
print('columns:', df.columns.tolist())
cols=[c for c in ['cell_id','centroid_x','centroid_y','cell_type','region_label'] if c in df.columns]
print(df[cols].head().to_string(index=False))
print('n=', len(df))
print('CA1 cells=', int((df['region_label']=='CA1').sum()))
print('Reactive Astrocyte in CA1=', int(((df['region_label']=='CA1') & (df['cell_type']=='Reactive Astrocyte')).sum()))
print('Pyramidal in CA1=', int(((df['region_label']=='CA1') & (df['cell_type']=='Pyramidal Neuron')).sum()))

print('\n--- region polygon names sample ---')
polys=load_region_polygons(slide)
print(sorted(polys.keys())[:12], '... total', len(polys))
PY
