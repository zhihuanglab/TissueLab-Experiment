set -e
python - <<'PY'
import sys
sys.path.insert(0, '/shared/lib')
from pathlib import Path
from shared_analysis.sea_ad_lfb import build_cell_table
slide = Path('/data/H19.33.004-A12-LFB.svs.zarr')
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print(cells.head(3).to_string())
print('\nUnique regions:', sorted([r for r in cells["region"].dropna().unique()])[:20])
ca1 = cells[cells['region']=='CA1']
print('\nCA1 n=', len(ca1))
print(ca1['cell_type'].value_counts().head(10).to_string())
PY
