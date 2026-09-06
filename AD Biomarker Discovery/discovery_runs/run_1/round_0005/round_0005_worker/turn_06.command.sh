set -e
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
from shared_analysis import build_cell_table, load_training_cohort
import shared_analysis.stats as stats
print('build_cell_table:', inspect.signature(build_cell_table))
print('load_training_cohort:', inspect.signature(load_training_cohort))
print('\nsource build_cell_table:')
print(inspect.getsource(build_cell_table)[:5000])
print('\nsource partial_correlation:')
print(inspect.getsource(stats.partial_correlation)[:2500])
print('\nsource leave_one_out_summary:')
print(inspect.getsource(stats.leave_one_out_summary)[:4000])
PY
echo '--- inspect one slide via build_cell_table columns/labels/regions ---'
python - <<'PY'
import sys
sys.path.insert(0, '/shared/lib')
from pathlib import Path
from shared_analysis import build_cell_table, load_training_cohort
cohort = load_training_cohort('/data')
row = cohort.iloc[0]
cells = build_cell_table(Path('/data')/row['slide_name'], include_regions=True, include_geometry=False)
print('slide', row['slide_name'], 'shape', cells.shape)
print('columns', cells.columns.tolist())
print('\nclass counts top:')
print(cells['class_label'].value_counts().head(10).to_string())
print('\nregion counts:')
print(cells['region_name'].value_counts(dropna=False).head(20).to_string())
print('\nhead:')
print(cells[['x_centroid','y_centroid','class_label','region_name']].head().to_string(index=False))
PY
