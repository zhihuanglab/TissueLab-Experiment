set -euo pipefail
python - <<'PY'
import sys, warnings, math
sys.path.insert(0,'/shared/lib')
from pathlib import Path
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from shared_analysis import build_cell_table, load_training_cohort

RUMS=[30.0,20.0]
MPP=0.503 # approx
cohort=load_training_cohort('/data')
rows=[]
for _,row in cohort.iterrows():
    donor=row['donor_id']; slide=Path('/data')/row['slide_name']
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cells=build_cell_table(slide, include_regions=True, include_geometry=False)
    ca1=cells[cells['region']=='CA1'][['x','y','cell_type']]
    ra=ca1[ca1['cell_type']=='Reactive Astrocyte'][['x','y']].to_numpy(float, copy=False)
    py=ca1[ca1['cell_type']=='Pyramidal Neuron'][['x','y']].to_numpy(float, copy=False)
    ly=ca1[ca1['cell_type']=='Lymphocyte'][['x','y']].to_numpy(float, copy=False)
    d={'donor_id':donor,'n_ca1':len(ca1),'n_ra':len(ra),'n_py':len(py),'n_ly':len(ly)}
    if len(ra)>0 and len(py)>0:
        ptree=cKDTree(py)
        dists,_=ptree.query(ra,k=1)
        peri = dists <= (30.0/MPP)
        d['n_peri_ra']=int(peri.sum())
        peri_xy=ra[peri]
        if len(peri_xy)>0 and len(ly)>0:
            ltree=cKDTree(ly)
            for rum in RUMS:
                ds,_=ltree.query(peri_xy,k=1)
                d[f'frac_ly_{int(rum)}um']=float((ds <= (rum/MPP)).mean())
        elif len(peri_xy)>0:
            for rum in RUMS:
                d[f'frac_ly_{int(rum)}um']=0.0
        else:
            for rum in RUMS:
                d[f'frac_ly_{int(rum)}um']=math.nan
    rows.append(d)
df=pd.DataFrame(rows)
print(df[['donor_id','n_ra','n_py','n_ly','n_peri_ra','frac_ly_30um','frac_ly_20um']].to_string(index=False))
print('\ncoverage:', df[['frac_ly_30um','frac_ly_20um']].notna().sum().to_dict())
print('nonzero lymphocyte donors:', int((df['n_ly']>0).sum()))
PY
