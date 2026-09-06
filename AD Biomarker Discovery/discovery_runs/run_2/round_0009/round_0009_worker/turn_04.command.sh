set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
python - <<'PY'
import sys
print('PYTHONPATH ok', sys.path[:5])
from shared_analysis import load_training_cohort, build_cell_table
from pathlib import Path
cohort = load_training_cohort('/data')
row = cohort.iloc[0]
slide_path = Path('/data') / row['slide_name']
print('sample slide:', row['donor_id'], slide_path)
cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
print('shape', cells.shape)
print('columns', list(cells.columns))
for col in ['cell_type','predicted_label','classification','region_name','region','x_centroid','y_centroid','centroid_x','centroid_y']:
    if col in cells.columns:
        print('\nCOL', col)
        print(cells[col].head().to_list())
label_col = next((c for c in ['cell_type','predicted_label','classification','label_name','class_label'] if c in cells.columns), None)
region_col = next((c for c in ['region_name','region','region_label'] if c in cells.columns), None)
print('\nlabel_col', label_col)
if label_col:
    print(cells[label_col].value_counts(dropna=False).head(20))
print('\nregion_col', region_col)
if region_col:
    print(cells[region_col].value_counts(dropna=False).head(20))
PY
