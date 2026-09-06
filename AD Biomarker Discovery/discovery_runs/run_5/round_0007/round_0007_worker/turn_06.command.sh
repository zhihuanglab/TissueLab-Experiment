set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0,'/shared/lib')
from pathlib import Path
import shared_analysis
print('shared_analysis ok:', shared_analysis.__file__)
from shared_analysis import load_training_cohort, build_cell_table
cohort = load_training_cohort('/data')
print('cohort shape:', cohort.shape)
row = cohort.iloc[0]
print('example donor/slide:', row['donor_id'], row['slide_name'])
slide_path = Path('/data')/row['slide_name']
cells = build_cell_table(slide_path, include_regions=True, include_geometry=True)
print('cell table shape:', cells.shape)
print('columns:', list(cells.columns))
for col in ['cell_type','region']:
    if col in cells.columns:
        print(f'\nvalue_counts for {col}:')
        print(cells[col].value_counts(dropna=False).head(20).to_string())
keep=[c for c in ['cell_type','region','centroid_x','centroid_y','x_centroid','y_centroid','area_um2','area'] if c in cells.columns]
print('\nhead:')
print(cells[keep].head().to_string())
PY
