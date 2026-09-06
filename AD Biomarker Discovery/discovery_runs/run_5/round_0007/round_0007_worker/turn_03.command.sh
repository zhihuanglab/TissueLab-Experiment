set -euo pipefail
python - <<'PY'
from pathlib import Path
import pandas as pd
from shared_analysis import load_training_cohort, build_cell_table
cohort = load_training_cohort('/data')
print('cohort columns:', list(cohort.columns))
row = cohort.iloc[0]
print('example donor/slide:', row['donor_id'], row['slide_name'])
slide_path = Path('/data')/row['slide_name']
cells = build_cell_table(slide_path, include_regions=True, include_geometry=True)
print('cell table shape:', cells.shape)
print('columns:', list(cells.columns))
for col in ['cell_type','region','x_centroid','y_centroid','area','label']:
    if col in cells.columns:
        print(f'[{col}] unique/sample:', cells[col].dropna().astype(str).head().tolist())
        if col in ('cell_type','region','label'):
            print(cells[col].value_counts(dropna=False).head(15))
PY
