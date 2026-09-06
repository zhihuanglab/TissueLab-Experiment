set -e
python - <<'PY'
import sys, pandas as pd, numpy as np
sys.path.insert(0,'/shared/lib')
from shared_analysis import load_training_cohort, build_cell_table
from scipy.spatial import cKDTree

data_root='/data'
cohort=load_training_cohort(data_root)
slide=cohort.iloc[0]['slide_name']
donor=cohort.iloc[0]['donor_id']
print('sample donor', donor, slide)
cells=build_cell_table(f'/data/{slide}', include_regions=True, include_geometry=False)
ca1=cells[cells['region']=='CA1']
print('total cells', len(cells), 'CA1 cells', len(ca1))
print(ca1['cell_type'].value_counts().to_string())
pyr=ca1[ca1['cell_type']=='Pyramidal Neuron'][['x','y']].to_numpy()
ra=ca1[ca1['cell_type']=='Reactive Astrocyte'][['x','y']].to_numpy()
print('pyramidal', len(pyr), 'reactive', len(ra))
if len(pyr) and len(ra):
    d,_=cKDTree(ra).query(pyr, k=1)
    for r in [4,6]:
        print('fraction within', r, float((d<=r).mean()))
    print('nearest distance quantiles', np.quantile(d,[0,.25,.5,.75,.9,.99]).tolist())
PY
