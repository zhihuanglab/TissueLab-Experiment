set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
from pathlib import Path
from shared_analysis import load_training_cohort, build_cell_table

print('build_cell_table source file:', inspect.getsourcefile(build_cell_table))
src = inspect.getsource(build_cell_table)
print('build_cell_table source head:')
print('\n'.join(src.splitlines()[:90]))

cohort = load_training_cohort('/data')
slide = Path('/data') / cohort.iloc[0]['slide_name']
print('\nTesting slide:', slide.name)
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print('shape', cells.shape)
print('columns', list(cells.columns))
for c in ['classification','cell_type','region','x_centroid','y_centroid','centroid_x','centroid_y']:
    if c in cells.columns:
        print('\nCOLUMN', c)
        print('nunique=', cells[c].nunique(dropna=False))
        print(cells[c].dropna().astype(str).value_counts().head(15))
PY
