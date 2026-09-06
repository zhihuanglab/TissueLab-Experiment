set -euo pipefail
python - <<'PY'
import sys, os
sys.path.insert(0, '/shared/lib')
from pathlib import Path
from shared_analysis import build_cell_table, load_region_polygons, load_training_cohort
data_root=Path('/data')
cohort=load_training_cohort(data_root)
slide_path=data_root/str(cohort.iloc[0]['slide_name'])
cells=build_cell_table(slide_path, include_regions=True, include_geometry=False)
print("CELL TABLE COLUMNS:", cells.columns.tolist())
print(cells.head().to_string())
print("region counts top:")
print(cells['region_label'].value_counts(dropna=False).head(15).to_string())
print("class counts top:")
print(cells['class_label'].value_counts(dropna=False).head(15).to_string())
polys=load_region_polygons(slide_path)
print("\nPOLYGON KEYS:", list(polys)[:20], "n=", len(polys))
for k,v in polys.items():
    print("example region", k, "type", type(v).__name__, "len", len(v) if hasattr(v,'__len__') else 'na')
    if hasattr(v,'__len__') and len(v):
        first=v[0]
        print("first entry type:", type(first).__name__)
        try:
            print("first 3 points:", first[:3])
        except Exception as e:
            print("cannot slice first entry", e)
    break
PY
