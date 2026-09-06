set -euo pipefail
echo '--- shared_analysis __init__ ---'
sed -n '1,220p' /shared/lib/shared_analysis/__init__.py
echo '--- sea_ad_lfb.py head ---'
sed -n '1,280p' /shared/lib/shared_analysis/sea_ad_lfb.py
echo '--- stats.py head ---'
sed -n '1,260p' /shared/lib/shared_analysis/stats.py
python - <<'PY'
import sys
sys.path.append('/shared/lib')
from pathlib import Path
from shared_analysis.sea_ad_lfb import load_training_cohort, build_cell_table
data_root=Path('/data')
cohort=load_training_cohort(data_root)
row=cohort.iloc[0]
slide_path=data_root/row['slide_name']
cells=build_cell_table(slide_path, include_regions=True, include_geometry=False)
print('OK built cells', cells.shape)
print('columns', cells.columns.tolist())
for col in ['cell_type','region_name','region_label','x_centroid','y_centroid']:
    if col in cells.columns:
        print(col, 'head', cells[col].head().tolist())
print('cell_type counts')
print(cells['cell_type'].value_counts().head(10).to_string())
print('regions', cells['region_name'].value_counts(dropna=False).head(10).to_string() if 'region_name' in cells.columns else 'no region_name')
PY
