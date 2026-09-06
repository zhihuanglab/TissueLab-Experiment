set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0, '/shared/lib')
from shared_analysis import load_training_cohort, build_cell_table
from pathlib import Path
cohort = load_training_cohort('/data')
print('cohort columns', list(cohort.columns))
row = cohort.iloc[0]
slide_path = Path('/data') / row['slide_name']
cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
print('shape', cells.shape)
print('columns', list(cells.columns))
label_col = next((c for c in ['cell_type','predicted_label','classification','label_name','class_label'] if c in cells.columns), None)
region_col = next((c for c in ['region_name','region','region_label'] if c in cells.columns), None)
xcol = next((c for c in ['x_centroid','centroid_x','x'] if c in cells.columns), None)
ycol = next((c for c in ['y_centroid','centroid_y','y'] if c in cells.columns), None)
print('label_col', label_col, 'region_col', region_col, 'xcol', xcol, 'ycol', ycol)
print('label counts top 15:\n', cells[label_col].value_counts(dropna=False).head(15).to_string())
print('region counts:\n', cells[region_col].value_counts(dropna=False).head(12).to_string())
print('sample rows:\n', cells[[xcol,ycol,label_col,region_col]].head().to_string())
PY
echo '--- shared_analysis __init__.py ---'
sed -n '1,240p' /shared/lib/shared_analysis/__init__.py
echo '--- artifacts.py head ---'
sed -n '1,260p' /shared/lib/shared_analysis/artifacts.py
echo '--- stats.py head ---'
sed -n '1,280p' /shared/lib/shared_analysis/stats.py
