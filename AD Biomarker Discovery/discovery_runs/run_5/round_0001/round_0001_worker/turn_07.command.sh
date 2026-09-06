set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
from pathlib import Path
from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort, load_class_lookup
print('signatures', inspect.signature(load_training_cohort), inspect.signature(build_cell_table))
src = inspect.getsource(build_cell_table)
print(src[:2200])
data_root=Path('/data')
cohort=load_training_cohort(data_root)
slide = data_root / cohort.iloc[0]['slide_name']
print('sample slide', slide.name)
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print('columns', list(cells.columns))
print(cells[['cell_type','region']].head(8).to_string())
print('top cell types')
print(cells['cell_type'].dropna().astype(str).value_counts().head(15).to_string())
print('regions')
print(cells['region'].dropna().astype(str).value_counts().to_string())
print('class lookup head', list(load_class_lookup(slide).items())[:10])
PY
