set -euo pipefail
python - <<'PY'
import sys, os
from pathlib import Path
sys.path.append('/shared/lib')
print('shared exists:', Path('/shared/lib').exists())
import shared_analysis
print('shared_analysis path:', shared_analysis.__file__)
from shared_analysis import build_cell_table, load_training_cohort
import pandas as pd

data_root=Path('/data')
cohort=load_training_cohort(data_root)
row=cohort.iloc[0]
slide=data_root/row['slide_name']
print('sample slide:', slide.name)

cells=build_cell_table(slide, include_regions=True, include_geometry=False)
print('shape', cells.shape)
print('columns', cells.columns.tolist())
print('cell labels', cells['cell_label'].value_counts().head(10).to_dict())
print('region labels', cells['region_label'].value_counts(dropna=False).head(12).to_dict())
ca1 = cells[cells['region_label']=='CA1']
print('CA1 n=', len(ca1))
print('CA1 labels', ca1['cell_label'].value_counts().head(10).to_dict())
PY
