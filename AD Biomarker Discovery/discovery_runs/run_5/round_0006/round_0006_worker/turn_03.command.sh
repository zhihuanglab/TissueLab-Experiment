set -euo pipefail
python - <<'PY'
from pathlib import Path
import pandas as pd
from shared_analysis import load_training_cohort, build_cell_table
cohort = load_training_cohort('/data')
row = cohort.iloc[0]
slide = Path('/data') / row['slide_name']
print('slide', slide)
cells = build_cell_table(slide, include_regions=True, include_geometry=True)
print('columns:', list(cells.columns))
print('n=', len(cells))
for col in ['cell_type','predicted_label','label','region','x_centroid','y_centroid','centroid_x','centroid_y','area','contour_area_um2']:
    if col in cells.columns:
        print(col, 'nonnull', cells[col].notna().sum(), 'unique sample', cells[col].dropna().astype(str).unique()[:10])
print('region counts top:')
for col in ['region','anatomical_region','region_label']:
    if col in cells.columns:
        print(col, cells[col].value_counts(dropna=False).head(10))
print('cell type counts top:')
for col in ['cell_type','predicted_label','label']:
    if col in cells.columns:
        print(col, cells[col].value_counts(dropna=False).head(10))
PY
