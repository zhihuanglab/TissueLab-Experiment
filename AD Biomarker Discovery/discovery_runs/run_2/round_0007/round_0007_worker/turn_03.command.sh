set -euo pipefail
python - <<'PY'
from pathlib import Path
import pandas as pd
import zarr
from shared_analysis import build_cell_table, load_training_cohort

data_root=Path('/data')
cohort=load_training_cohort(data_root)
row=cohort.iloc[0]
slide_path=data_root/row['slide_name']
print('sample slide:', slide_path)

print('\n--- zarr tree (depth 2) ---')
root=zarr.open(str(slide_path), mode='r')
def walk(g, prefix='', depth=0, maxdepth=2):
    if depth>maxdepth: return
    for k,v in g.groups():
        print(prefix+k+'/')
        walk(v, prefix+'  ', depth+1, maxdepth)
    for k,v in g.arrays():
        try: shape=v.shape; dtype=v.dtype
        except Exception: shape='?'; dtype='?'
        print(prefix+f'{k} shape={shape} dtype={dtype}')
walk(root)

print('\n--- build_cell_table sample ---')
cells=build_cell_table(slide_path, include_regions=True, include_geometry=False)
print(cells.head().to_string())
print('\ncolumns:', list(cells.columns))
print('\nregion counts top:')
print(cells['region_name'].value_counts(dropna=False).head(12).to_string())
print('\nlabel counts:')
label_col=[c for c in cells.columns if 'label' in c.lower() or 'class' in c.lower()]
print(label_col)
for c in label_col:
    print(cells[c].value_counts().head(10).to_string())
PY
