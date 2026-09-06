set -euo pipefail
python - <<'PY'
import sys
sys.path.append('/shared/lib')
from pathlib import Path
from shared_analysis import build_cell_table, load_training_cohort
cohort = load_training_cohort('/data')
slide = Path('/data') / cohort.iloc[0]['slide_name']
print('slide:', slide.name)
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print('columns:', cells.columns.tolist())
print('n=', len(cells))
print('cell types:', cells['cell_type'].value_counts().head(20).to_dict())
print('regions:', cells['region'].value_counts(dropna=False).to_dict())
print(cells.head().to_string())
PY
