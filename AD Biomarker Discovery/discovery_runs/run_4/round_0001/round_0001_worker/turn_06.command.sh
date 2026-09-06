set -euo pipefail
python - <<'PY'
import sys, math
sys.path.insert(0, '/shared/lib')
from pathlib import Path
import numpy as np, pandas as pd
from shared_analysis import build_cell_table, load_region_polygons, load_training_cohort
data_root=Path('/data')
cohort=load_training_cohort(data_root)
rows=[]
for _,r in cohort.iterrows():
    slide_path=data_root/str(r['slide_name'])
    cells=build_cell_table(slide_path, include_regions=True, include_geometry=False)
    ca1=cells[cells['region']=='CA1']
    ca1_total=len(ca1)
    ca1_pyr=int((ca1['cell_type']=='Pyramidal Neuron').sum())
    polys=load_region_polygons(slide_path)
    area=np.nan
    if 'CA1' in polys:
        area=sum(0.5*abs(np.dot(p[:,0], np.roll(p[:,1], -1)) - np.dot(p[:,1], np.roll(p[:,0], -1))) for p in polys['CA1'])
    rows.append((r['donor_id'], ca1_total, ca1_pyr, area))
df=pd.DataFrame(rows, columns=['donor_id','ca1_total','ca1_pyr','ca1_area'])
df['fraction']=df['ca1_pyr']/df['ca1_total']
df['density']=df['ca1_pyr']/df['ca1_area']
print(df.describe().to_string())
print("\nlowest ca1_total:")
print(df.sort_values('ca1_total').head(10).to_string(index=False))
print("\nlowest ca1_pyr:")
print(df.sort_values('ca1_pyr').head(10).to_string(index=False))
PY
