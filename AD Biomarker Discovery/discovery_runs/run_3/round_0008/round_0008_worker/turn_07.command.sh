set -euo pipefail
python - <<'PY'
import sys
sys.path.append('/shared/lib')
from pathlib import Path
from shared_analysis import load_training_cohort, build_cell_table
cohort=load_training_cohort('/data')
row=cohort.iloc[0]
print('sample donor', row['donor_id'], row['slide_name'])
cells=build_cell_table(Path('/data')/row['slide_name'], include_regions=True, include_geometry=False)
print(cells.head().to_string())
print('columns:', list(cells.columns))
for col in ['cell_type','predicted_class','region_label','region_name','x_centroid','y_centroid','centroid_x','centroid_y']:
    if col in cells.columns:
        vals=cells[col].dropna().astype(str).unique()[:20]
        print(col, vals)
print('n=', len(cells))
PY
