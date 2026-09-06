set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0,'/shared/lib')
from pathlib import Path
from shared_analysis import load_training_cohort, build_cell_table
cohort = load_training_cohort('/data')
row = cohort.iloc[0]
cells = build_cell_table(Path('/data') / row['slide_name'], include_regions=True, include_geometry=False)
print(sorted(set(cells['cell_type'].tolist())))
print(sorted(set([x for x in cells['region'].tolist() if x==x])))
PY
