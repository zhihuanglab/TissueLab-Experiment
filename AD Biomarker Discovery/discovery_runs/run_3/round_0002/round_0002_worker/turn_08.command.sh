set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0,'/shared/lib')
from pathlib import Path
import zarr
from shared_analysis import load_training_cohort, build_cell_table

data_root=Path('/data')
cohort=load_training_cohort(data_root)
row=cohort.iloc[0]
slide=row['slide_name']
print('slide', slide)
z = zarr.open(str(data_root/slide), mode='r')
print('top groups:', list(z.group_keys()))
for g in z.group_keys():
    grp=z[g]
    print(f'group {g}: subgroups={list(grp.group_keys())[:10]} arrays={list(grp.array_keys())[:10]}')
cells = build_cell_table(data_root/slide, include_regions=True, include_geometry=False)
print('\nCELL TABLE SHAPE:', cells.shape)
print('COLUMNS:', list(cells.columns))
for c in ['cell_type','predicted_label','classification','region','x_centroid','y_centroid']:
    if c in cells.columns:
        print(f'\n{c} value counts head:')
        print(cells[c].dropna().astype(str).value_counts().head(10).to_string())
print('\nHEAD:')
print(cells.head().to_string())
PY
