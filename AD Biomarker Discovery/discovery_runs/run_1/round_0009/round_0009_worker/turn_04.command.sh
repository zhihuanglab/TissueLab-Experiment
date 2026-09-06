set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
python - <<'PY'
from pathlib import Path
from shared_analysis import build_cell_table, load_training_cohort

cohort = load_training_cohort('/data')
row = cohort.iloc[0]
slide = Path('/data')/row['slide_name']
print('slide', slide)
cells = build_cell_table(slide, include_regions=True, include_geometry=True)
print('shape', cells.shape)
print('columns', list(cells.columns))
print('cell_type counts head', cells['cell_type'].value_counts().head(12).to_dict())
print('region counts', cells['region_label'].value_counts().to_dict())
print(cells[['cell_type','region_label','centroid_x','centroid_y']].head().to_string())
geom_cols = [c for c in cells.columns if 'contour' in c.lower() or 'area' in c.lower()]
print('geom cols', geom_cols)
for c in geom_cols:
    v = cells[c].iloc[0]
    print(c, type(v), getattr(v,'shape',None), (len(v) if hasattr(v,'__len__') and not isinstance(v,(str,bytes)) else None))
PY
