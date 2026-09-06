set -euo pipefail
python - <<'PY'
from pathlib import Path
import inspect, json
from shared_analysis import build_cell_table, load_training_cohort
import pandas as pd
import zarr

print('build_cell_table signature:', inspect.signature(build_cell_table))
cohort=load_training_cohort('/data')
slide=Path('/data')/cohort.iloc[0]['slide_name']
cells=build_cell_table(slide, include_regions=True, include_geometry=True)
print('cell_table columns:', list(cells.columns))
print(cells.head(3).to_dict(orient='records'))
print('n=', len(cells), 'slide=', slide.name)

print('\nzarr tree-ish top level:')
root=zarr.open(str(slide), mode='r')
print(list(root.group_keys()))
for g in root.group_keys():
    grp=root[g]
    print('group', g, 'arrays', list(grp.array_keys())[:10], 'subgroups', list(grp.group_keys())[:10])
PY
