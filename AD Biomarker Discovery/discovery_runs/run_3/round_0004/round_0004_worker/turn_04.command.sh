set -euo pipefail
python - <<'PY'
import sys
sys.path.append('/shared/lib')
from pathlib import Path
from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort

data_root = Path('/data')
cohort = load_training_cohort(data_root)
slide_name = cohort.iloc[0]['slide_name']
print('slide_name', slide_name)
cells = build_cell_table(data_root/slide_name, include_regions=True, include_geometry=False)
print('columns:', cells.columns.tolist())
print('shape:', cells.shape)
for col in ['class_label','region_label','centroid_x','centroid_y']:
    if col in cells.columns:
        print(col, cells[col].head().tolist())
print('class counts:', cells['class_label'].value_counts().head(10).to_dict())
print('region counts:', cells['region_label'].value_counts(dropna=False).head(10).to_dict())
PY
