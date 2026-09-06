set -euo pipefail
python - <<'PY'
import inspect, textwrap, zarr, os, json
from pathlib import Path
from shared_analysis.sea_ad_lfb import build_cell_table

print('--- build_cell_table source snippet ---')
src = inspect.getsource(build_cell_table)
print(src[:5000])

slide = Path('/data/H19.33.004-A12-LFB.svs.zarr')
print('\n--- sample zarr tree top-level ---')
root = zarr.open(slide, mode='r')
print(root.tree())

print('\n--- sample cell table columns/head ---')
df = build_cell_table(slide, include_regions=True, include_geometry=False)
print(df.columns.tolist())
print(df[['cell_id','centroid_x','centroid_y','cell_type','region_label']].head().to_string(index=False))
print('n=', len(df), 'regions=', df['region_label'].value_counts(dropna=False).head(10).to_dict())
PY
