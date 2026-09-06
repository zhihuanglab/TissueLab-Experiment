set -euo pipefail
python - <<'PY'
from pathlib import Path
import inspect
from shared_analysis import build_cell_table, load_training_cohort
print('build_cell_table source:')
print(inspect.getsource(build_cell_table))
cohort = load_training_cohort('/data')
slide = Path('/data') / cohort.iloc[0]['slide_name']
print('\nslide:', slide)
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print('columns:', cells.columns.tolist())
print('n=', len(cells))
for col in ['cell_type','region','centroid_x','centroid_y']:
    if col in cells.columns:
        print(col, 'sample:', cells[col].head().tolist())
print('unique cell_type sample:', cells['cell_type'].dropna().value_counts().head(15).to_dict())
print('unique region:', cells['region'].dropna().value_counts().to_dict())
PY
