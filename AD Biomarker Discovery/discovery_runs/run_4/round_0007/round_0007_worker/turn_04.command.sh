set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
python - <<'PY'
from pathlib import Path
import pandas as pd
from shared_analysis import load_training_cohort, build_cell_table
data_root=Path('/data')
cohort=load_training_cohort(data_root)
print('cohort columns:', list(cohort.columns))
row=cohort.iloc[0]
slide_path=data_root/row['slide_name']
print('slide:', slide_path.name)
cells=build_cell_table(slide_path, include_regions=True, include_geometry=True)
print('columns:', list(cells.columns))
print('n_cells:', len(cells))
print('cell_type counts top 12:')
print(cells['cell_type'].value_counts().head(12).to_string())
print('region counts:')
print(cells['region_label'].dropna().value_counts().to_string())
print('sample rows:')
cols=[c for c in ['cell_type','x_centroid','y_centroid','region_label','contour','contour_x','contour_y'] if c in cells.columns]
print(cells[cols].head(3).to_string())
geom_cols=[c for c in cells.columns if 'contour' in c.lower() or 'area' in c.lower()]
print('geom-like cols:', geom_cols)
PY
