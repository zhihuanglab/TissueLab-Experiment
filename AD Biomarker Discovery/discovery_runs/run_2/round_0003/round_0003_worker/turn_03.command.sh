set -euo pipefail
python - <<'PY'
from pathlib import Path
import pandas as pd
from shared_analysis import build_cell_table, load_training_cohort

data_root=Path('/data')
cohort=load_training_cohort(data_root)
print('cohort columns:', cohort.columns.tolist())
row=cohort.iloc[0]
slide=data_root/row['slide_name']
print('sample slide:', slide)

cells=build_cell_table(slide, include_regions=True, include_geometry=False)
print('cell table shape:', cells.shape)
print('cell table columns:', cells.columns.tolist())
print('cell_label counts head:')
print(cells['cell_label'].value_counts().head(10))
print('region counts head:')
print(cells['region_label'].value_counts(dropna=False).head(15))
print('CA1 labels cross-tab head:')
print(pd.crosstab(cells['cell_label'], cells['region_label']).get('CA1', pd.Series(dtype=int)).sort_values(ascending=False).head(10))
print('null region frac:', cells['region_label'].isna().mean())
PY
