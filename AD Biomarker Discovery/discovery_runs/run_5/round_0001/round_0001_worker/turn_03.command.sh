set -euo pipefail
python - <<'PY'
from pathlib import Path
from shared_analysis import load_training_cohort, build_cell_table
import pandas as pd

data_root=Path('/data')
cohort=load_training_cohort(data_root)
slide = data_root / cohort.iloc[0]['slide_name']
print('sample slide:', slide)
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print('columns:', list(cells.columns))
for col in ['cell_type','region','donor_id']:
    if col in cells.columns:
        vals = cells[col].dropna().astype(str)
        print(f'{col} unique sample:', vals.unique()[:20], 'n_unique=', vals.nunique())
print('nrows', len(cells))
print(cells[['cell_type','region']].head(10).to_string())
PY
