set -e
python - <<'PY'
import sys
sys.path.insert(0,'/shared/lib')
from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort
from pathlib import Path
import pandas as pd

cohort = load_training_cohort('/data')
print('cohort columns:', list(cohort.columns))
slide = Path('/data') / cohort.iloc[0]['slide_name']
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print('cells shape:', cells.shape)
print('cell types:', sorted(cells['cell_type'].dropna().unique().tolist()))
print('regions:', sorted(cells['region'].dropna().unique().tolist()))
print('CA1 counts by cell_type head:')
print(cells.loc[cells['region']=='CA1','cell_type'].value_counts().head(10))
print('fraction unassigned region:', cells['region'].isna().mean())
PY
