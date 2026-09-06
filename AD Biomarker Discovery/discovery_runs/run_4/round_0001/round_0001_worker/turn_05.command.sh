set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0, '/shared/lib')
from pathlib import Path
from shared_analysis import build_cell_table, load_region_polygons, load_training_cohort
data_root=Path('/data')
cohort=load_training_cohort(data_root)
slide_path=data_root/str(cohort.iloc[0]['slide_name'])
cells=build_cell_table(slide_path, include_regions=True, include_geometry=False)
print("non-null region counts:")
print(cells['region'].value_counts(dropna=False).head(12).to_string())
ca1 = cells[cells['region']=='CA1']
print("\nCA1 rows:", len(ca1), "pyramidal:", int((ca1['cell_type']=='Pyramidal Neuron').sum()))
polys=load_region_polygons(slide_path)
print("\nPOLYGON REGION KEYS:", list(polys.keys()))
for reg in ['CA1','CA2','DG']:
    v=polys.get(reg)
    print("\nregion", reg, "type", type(v).__name__)
    if v is None:
        continue
    print("len", len(v))
    first=v[0]
    print("first poly type", type(first).__name__, "len", len(first))
    print("first 5 pts", first[:5])
PY
