set -euo pipefail
python - <<'PY'
from pathlib import Path
import pandas as pd
from shared_analysis import build_cell_table, load_training_cohort

cohort = load_training_cohort('/data')
row = cohort.iloc[0]
slide = Path('/data')/row['slide_name']
print('slide', slide)
cells = build_cell_table(slide, include_regions=True, include_geometry=True)
print('shape', cells.shape)
print('columns', list(cells.columns))
for col in ['cell_type','region_label','centroid_x','centroid_y','contour_x','contour_y','area_px2']:
    print(col, 'present=', col in cells.columns, 'dtype=', cells[col].dtype if col in cells.columns else None)
print('unique cell types sample', cells['cell_type'].dropna().value_counts().head(15).to_dict())
print('unique regions', cells['region_label'].dropna().value_counts().to_dict())
print(cells[['cell_type','region_label','centroid_x','centroid_y']].head().to_string())
# Show example geometry-related entries
geom_cols = [c for c in cells.columns if 'contour' in c.lower() or 'area' in c.lower()]
print('geom cols', geom_cols)
if 'contour_x' in cells.columns:
    print('first contour_x type', type(cells['contour_x'].iloc[0]), 'len', len(cells['contour_x'].iloc[0]))
    print('first contour_y type', type(cells['contour_y'].iloc[0]), 'len', len(cells['contour_y'].iloc[0]))
PY
