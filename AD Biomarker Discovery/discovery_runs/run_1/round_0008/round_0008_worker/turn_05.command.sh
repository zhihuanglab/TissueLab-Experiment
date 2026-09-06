set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
python - <<'PY'
from pathlib import Path
from shared_analysis import build_cell_table, load_training_cohort

data_root=Path('/data')
cohort=load_training_cohort(data_root)
row=cohort.iloc[0]
slide=data_root/row['slide_name']
print('slide', slide.name)
cells=build_cell_table(slide, include_regions=True, include_geometry=True)
print('shape', cells.shape)
print('columns', list(cells.columns))
for col in ['cell_type','predicted_label','label','region','x_centroid','y_centroid','contour_area','area']:
    if col in cells.columns:
        s=cells[col]
        print('\n', col, 'nonnull', int(s.notna().sum()))
        print(' sample:', s.dropna().astype(str).head().tolist())
if 'cell_type' in cells.columns:
    print('\ncell_type counts:', cells['cell_type'].value_counts().head(12).to_dict())
if 'region' in cells.columns:
    print('\nregion counts:', cells['region'].value_counts(dropna=False).head(12).to_dict())
PY
