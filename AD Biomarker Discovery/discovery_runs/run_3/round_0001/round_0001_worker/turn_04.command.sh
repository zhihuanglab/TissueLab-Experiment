set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
python - <<'PY'
from pathlib import Path
from shared_analysis.sea_ad_lfb import load_training_cohort, build_cell_table
import pandas as pd, json
data_root=Path('/data')
cohort=load_training_cohort(data_root)
print('cohort shape', cohort.shape)
row=cohort.iloc[0]
print('example donor', row['donor_id'], row['slide_name'])
slide_path=data_root/row['slide_name']
print('slide exists', slide_path.exists())
# inspect zarr tree shallow
import zarr
root=zarr.open(str(slide_path), mode='r')
print('root groups:', list(root.group_keys()))
for g in root.group_keys():
    grp=root[g]
    print(' group', g, 'arrays', list(grp.array_keys())[:10], 'groups', list(grp.group_keys())[:10])

cells=build_cell_table(slide_path, include_regions=True, include_geometry=False)
print('cell table columns', cells.columns.tolist())
for col in ['classification_label','predicted_label','region_label','region']:
    if col in cells.columns:
        print('\\nvalue counts for', col)
        print(cells[col].value_counts(dropna=False).head(20).to_string())
PY
