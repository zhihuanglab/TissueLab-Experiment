set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
python - <<'PY'
import inspect, zarr
from pathlib import Path
from shared_analysis.sea_ad_lfb import build_cell_table

print('--- build_cell_table source snippet ---')
src = inspect.getsource(build_cell_table)
print(src[:4500])

slide = Path('/data/H19.33.004-A12-LFB.svs.zarr')
print('\n--- sample zarr groups ---')
root = zarr.open(slide, mode='r')
print(root.tree()[:3000] if hasattr(root.tree(),'__str__') else root.tree())

print('\n--- sample cell table columns/head ---')
df = build_cell_table(slide, include_regions=True, include_geometry=False)
cols=[c for c in ['cell_id','centroid_x','centroid_y','cell_type','region_label','classification_label','class_label'] if c in df.columns]
print(df.columns.tolist())
print(df[cols].head().to_string(index=False))
print('n=', len(df))
print('region counts top:', df['region_label'].value_counts(dropna=False).head(10).to_dict())
print('cell_type counts top:', df['cell_type'].value_counts(dropna=False).head(10).to_dict())
PY
