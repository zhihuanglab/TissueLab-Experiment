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
cells=build_cell_table(slide_path, include_regions=True, include_geometry=True)
print('columns:', list(cells.columns))
print('n_cells:', len(cells))
print('cell_type value counts head:')
print(cells['cell_type'].value_counts().head(15))
print('regions:', cells['region_label'].dropna().value_counts().head(10))
print('sample rows:')
print(cells[['cell_type','x_centroid','y_centroid','region_label']].head().to_string())
geom_cols=[c for c in cells.columns if 'contour' in c.lower() or 'area' in c.lower()]
print('geom-like cols:', geom_cols)
if 'contour' in cells.columns:
    print('contour sample type:', type(cells['contour'].iloc[0]), cells['contour'].iloc[0][:3])
PY
