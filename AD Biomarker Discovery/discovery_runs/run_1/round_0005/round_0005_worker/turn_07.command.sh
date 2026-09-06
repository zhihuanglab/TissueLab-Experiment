set -e
python - <<'PY'
import sys, warnings
sys.path.insert(0, '/shared/lib')
warnings.filterwarnings('ignore')
from pathlib import Path
from shared_analysis import build_cell_table, load_training_cohort
cohort = load_training_cohort('/data')
row = cohort.iloc[0]
cells = build_cell_table(Path('/data')/row['slide_name'], include_regions=True, include_geometry=False)
print('slide', row['slide_name'], 'shape', cells.shape)
print('class counts top:')
print(cells['cell_type'].value_counts().to_string())
print('\nregion counts:')
print(cells['region'].value_counts(dropna=False).to_string())
print('\nhead:')
print(cells[['x','y','cell_type','region']].head().to_string(index=False))
PY
