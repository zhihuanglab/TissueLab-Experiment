set -euo pipefail
python - <<'PY'
from shared_analysis import load_training_cohort, build_cell_table
from pathlib import Path
cohort = load_training_cohort('/data')
row = cohort.iloc[0]
slide_path = Path('/data') / row['slide_name']
print('sample slide:', row['donor_id'], slide_path)
cells = build_cell_table(slide_path, include_regions=True, include_geometry=False)
print('shape', cells.shape)
print('columns', list(cells.columns))
for col in ['cell_type','predicted_label','classification','region_name','region','x_centroid','y_centroid','centroid_x','centroid_y']:
    if col in cells.columns:
        print('\nCOL', col)
        print(cells[col].head().to_list())
print('\nregion unique sample:', cells[[c for c in cells.columns if "region" in c.lower()][:4]].head())
label_col = None
for cand in ['cell_type','predicted_label','classification','label_name','class_label']:
    if cand in cells.columns:
        label_col = cand
        break
print('\nlabel_col', label_col)
if label_col:
    vc = cells[label_col].value_counts(dropna=False).head(20)
    print(vc)
region_col = None
for cand in ['region_name','region','region_label']:
    if cand in cells.columns:
        region_col = cand
        break
print('\nregion_col', region_col)
if region_col:
    print(cells[region_col].value_counts(dropna=False).head(20))
PY
