set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0,'/shared/lib')
from pathlib import Path
import inspect
from shared_analysis import build_cell_table, load_training_cohort
import zarr

print('build_cell_table signature:', inspect.signature(build_cell_table))
cohort=load_training_cohort('/data')
slide=Path('/data')/cohort.iloc[0]['slide_name']
cells=build_cell_table(slide, include_regions=True, include_geometry=True)
print('cell_table columns:', list(cells.columns))
print('head row subset:', {k:v for k,v in cells.head(1).to_dict(orient='records')[0].items() if k in ['donor_id','slide_name','cell_index','cell_type','region_label','centroid_x','centroid_y','contour_area_px']})
print('n=', len(cells), 'slide=', slide.name)

print('\nzarr groups:')
root=zarr.open(str(slide), mode='r')
for g in root.group_keys():
    grp=root[g]
    print('group', g, 'arrays', list(grp.array_keys())[:10], 'subgroups', list(grp.group_keys())[:10])
PY
