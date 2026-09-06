set -euo pipefail
python - <<'PY'
from pathlib import Path
import pandas as pd
from shared_analysis import load_training_cohort, build_cell_table
data_root=Path('/data')
cohort=load_training_cohort(data_root)
row=cohort.iloc[0]
slide_path=data_root/row['slide_name']
print('slide:', slide_path)
cells=build_cell_table(slide_path, include_regions=True, include_geometry=False)
print('shape', cells.shape)
print('columns', list(cells.columns))
for col in [c for c in cells.columns if 'label' in c.lower() or 'class' in c.lower() or 'region' in c.lower()]:
    vals=cells[col].dropna().astype(str)
    uniq=vals.unique()[:20]
    print('\nCOL',col,'nuniq',vals.nunique(),'sample',uniq)
print('\nregion counts top:')
if 'region_label' in cells.columns:
    print(cells['region_label'].value_counts(dropna=False).head(20))
elif 'region' in cells.columns:
    print(cells['region'].value_counts(dropna=False).head(20))
print('\ncell type counts top:')
for cand in ['class_label','cell_type','label','predicted_label','classification_label']:
    if cand in cells.columns:
        print(cand, cells[cand].value_counts().head(20))
PY
