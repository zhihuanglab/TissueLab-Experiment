set -euo pipefail
python - <<'PY'
from pathlib import Path
import inspect
from shared_analysis import build_cell_table, load_training_cohort
import shared_analysis.sea_ad_lfb as m

print('build_cell_table source:')
print(inspect.getsource(build_cell_table))
print('\nload_region_polygons source:')
print(inspect.getsource(m.load_region_polygons))

cohort = load_training_cohort('/data')
slide = Path('/data') / cohort.iloc[0]['slide_name']
print('\nSample slide:', slide)
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print('columns:', list(cells.columns))
for col in ['region','cell_type','classification_label','predicted_label']:
    if col in cells.columns:
        print(col, 'nunique=', cells[col].nunique(dropna=False))
        print(cells[col].value_counts(dropna=False).head(15).to_string())
PY
