set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
python - <<'PY'
from pathlib import Path
import pandas as pd
from shared_analysis.sea_ad_lfb import load_training_cohort, build_cell_table
data_root=Path('/data')
cohort=load_training_cohort(data_root)
row=cohort.iloc[0]
slide_path=data_root/row['slide_name']
print('slide_path', slide_path)
cells=build_cell_table(slide_path, include_regions=True, include_geometry=False)
print('columns:', list(cells.columns))
print('n_cells:', len(cells))
for col in ['cell_type','region_name','x_centroid','y_centroid','label']:
    if col in cells.columns:
        print('\n', col, cells[col].head().tolist())
if 'cell_type' in cells.columns:
    print('\ncell types head counts:')
    print(cells['cell_type'].value_counts().head(15).to_string())
if 'region_name' in cells.columns:
    print('\nregions counts:')
    print(cells['region_name'].value_counts(dropna=False).head(12).to_string())
PY
