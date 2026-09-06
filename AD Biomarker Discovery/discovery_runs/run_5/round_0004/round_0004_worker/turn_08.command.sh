set -euo pipefail
echo '--- sea_ad_lfb helper source ---'
sed -n '1,260p' /shared/lib/shared_analysis/sea_ad_lfb.py
echo '--- stats helper source ---'
sed -n '1,260p' /shared/lib/shared_analysis/stats.py
echo '--- one-slide data sanity via direct sys.path append ---'
python - <<'PY'
import sys
sys.path.append('/shared/lib')
from pathlib import Path
import pandas as pd
from shared_analysis.sea_ad_lfb import build_cell_table, load_region_polygons
cohort = pd.read_csv('/data/training_cohort.csv')
slide = Path('/data') / cohort.iloc[0]['slide_name']
print('slide', slide)
polys = load_region_polygons(slide)
print('region keys', list(polys.keys())[:10], 'n=', len(polys))
cells = build_cell_table(slide, include_regions=True, include_geometry=True)
print('n cells', len(cells))
print(cells.columns.tolist())
show=[c for c in ['cell_label','region_label','centroid_x','centroid_y','contour_area'] if c in cells.columns]
print(cells[show].head().to_string())
print('regions value counts head')
print(cells['region_label'].value_counts(dropna=False).head(12).to_string())
print('dtype cell_label', cells['cell_label'].dtype)
PY
