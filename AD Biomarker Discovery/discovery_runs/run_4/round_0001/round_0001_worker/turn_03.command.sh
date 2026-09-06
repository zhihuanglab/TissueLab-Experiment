set -euo pipefail
python - <<'PY'
from pathlib import Path
import pandas as pd
from shared_analysis import build_cell_table, load_region_polygons, load_training_cohort
data_root=Path('/data')
cohort=load_training_cohort(data_root)
slide_path=data_root/str(cohort.iloc[0]['slide_name'])
cells=build_cell_table(slide_path, include_regions=True, include_geometry=False)
print("CELL TABLE COLUMNS:", cells.columns.tolist())
print(cells.head().to_string())
print("region counts top:", cells['region_label'].value_counts(dropna=False).head(15).to_string())
print("class counts top:", cells['class_label'].value_counts(dropna=False).head(15).to_string())
polys=load_region_polygons(slide_path)
print("\nPOLYGON KEYS:", list(polys)[:20], "n=", len(polys))
for k,v in polys.items():
    print("example region", k, "type", type(v).__name__, "len", len(v) if hasattr(v,'__len__') else 'na')
    if hasattr(v,'__len__') and len(v):
        print("first entry type:", type(v[0]).__name__)
        print("first 3 points:", v[0][:3] if hasattr(v[0],'__getitem__') else v[0])
    break
PY
