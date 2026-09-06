set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0,'/shared/lib')
from pathlib import Path
from shared_analysis import load_training_cohort, build_cell_table
cohort = load_training_cohort('/data')
row = cohort.iloc[0]
slide = Path('/data') / row['slide_name']
print('slide', slide.name)
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print('columns:', list(cells.columns))
print('n=', len(cells))
for col in ['cell_type','predicted_label','label','region','anatomical_region','region_label','x_centroid','y_centroid','centroid_x','centroid_y']:
    if col in cells.columns:
        vals = cells[col]
        print('\nCOL', col)
        print('dtype', vals.dtype)
        print('head', vals.head().tolist())
        print('nunique', vals.nunique(dropna=True))
        if vals.dtype == object:
            print('top', vals.value_counts(dropna=False).head(10).to_dict())
PY
