set -euo pipefail
python - <<'PY'
import sys
sys.path.append('/shared/lib')
from pathlib import Path
from shared_analysis import load_training_cohort, build_cell_table
cohort=load_training_cohort('/data')
for i in [0,1]:
    row=cohort.iloc[i]
    cells=build_cell_table(Path('/data')/row['slide_name'], include_regions=True, include_geometry=False)
    print('\nslide', row['slide_name'])
    print(cells['region'].value_counts(dropna=False).head(15).to_string())
    sel=cells[cells['region']=='CA1']
    print('CA1 cell types:\n', sel['cell_type'].value_counts().head(10).to_string())
PY
