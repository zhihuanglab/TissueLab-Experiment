set -euo pipefail
python - <<'PY'
from pathlib import Path
import pandas as pd
from shared_analysis import load_training_cohort, build_cell_table
data_root=Path('/data')
cohort=load_training_cohort(data_root)
print('cohort shape', cohort.shape)
row=cohort.iloc[0]
print('example donor', row['donor_id'], row['slide_name'])
slide_path=data_root/row['slide_name']
print('slide exists', slide_path.exists())
# inspect zarr tree shallow
try:
    import zarr
    root=zarr.open(str(slide_path), mode='r')
    print('zarr groups:', list(root.group_keys()))
    for g in root.group_keys():
        grp=root[g]
        print(' group', g, 'arrays', list(grp.array_keys())[:10], 'groups', list(grp.group_keys())[:10])
except Exception as e:
    print('zarr inspect failed', e)

cells=build_cell_table(slide_path, include_regions=True, include_geometry=False)
print('cell table columns', cells.columns.tolist())
print('n cells', len(cells))
for col in ['classification_label','predicted_label','region_label','region']:
    if col in cells.columns:
        print('value counts for', col)
        print(cells[col].value_counts(dropna=False).head(20))
PY
