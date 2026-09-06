set -euo pipefail
echo '--- /shared/lib ---'
find /shared/lib -maxdepth 2 -type f | sed -n '1,120p'
echo '--- inspect shared_analysis exports and sample cell table ---'
python - <<'PY'
import sys
sys.path.insert(0, '/shared/lib')
from pathlib import Path
import shared_analysis
print('shared_analysis imported from', shared_analysis.__file__)
from shared_analysis import load_training_cohort, build_cell_table
data_root=Path('/data')
cohort=load_training_cohort(data_root)
row=cohort.iloc[0]
slide_path=data_root/row['slide_name']
cells=build_cell_table(slide_path, include_regions=True, include_geometry=False)
print('shape', cells.shape)
print('columns', list(cells.columns))
for col in ['region_name','region_label','region','classification','class_label','cell_type','label','x','y','centroid_x','centroid_y']:
    if col in cells.columns:
        vals=cells[col]
        print('\nCOL', col)
        if vals.dtype=='O':
            print(vals.value_counts(dropna=False).head(10))
        else:
            print(vals.head())
PY
