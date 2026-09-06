set -euo pipefail
python - <<'PY'
import sys, warnings
from pathlib import Path
sys.path.append('/shared/lib')
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort
from shared_analysis.stats import partial_correlation, residualized_loo_predictive_correlation
try:
    import openslide
except Exception:
    openslide = None
warnings.filterwarnings('ignore')
cohort=load_training_cohort('/data').copy()
cohort['sex_binary']=cohort['sex'].map({'Female':0.0,'Male':1.0})
rows=[]
for row in cohort.itertuples(index=False):
    zarr=Path('/data')/row.slide_name
    svs=zarr.with_suffix('')
    cells=build_cell_table(zarr, include_regions=True, include_geometry=False)
    ca1=cells['region'].eq('CA1')
    ct=cells['cell_type']
    pyr=cells.loc[ca1 & ct.eq('Pyramidal Neuron'), ['x','y']]
    lineage=cells.loc[ca1 & ct.isin(['Astrocyte','Reactive Astrocyte']), ['x','y','cell_type']]
    if openslide is not None:
        try:
            slide=openslide.OpenSlide(str(svs))
            mpp=float(slide.properties.get('openslide.mpp-x') or slide.properties.get('aperio.MPP') or 0.503)
        except Exception:
            mpp=0.503
    else:
        mpp=0.503
    out={'donor_id':row.donor_id,'outcome':row.slope_zmem0,'max_age_vis':row.max_age_vis,'braak_numeric':row.braak_numeric,'cerad_ordinal':row.cerad_ordinal,'sex_binary':row.sex_binary,
         'ca1_pyr':len(pyr),'ca1_lineage':len(lineage)}
    if len(pyr)==0 or len(lineage)==0:
        out.update({'f20':np.nan,'f30':np.nan,'n20':0,'n30':0,'r20':0,'r30':0,'a20':0,'a30':0})
    else:
        tree=cKDTree(pyr[['x','y']].to_numpy(float))
        dist,_=tree.query(lineage[['x','y']].to_numpy(float),k=1,workers=-1)
        for um in [20,30]:
            near=dist <= (um/mpp)
            sub=lineage.loc[near,'cell_type']
            r=(sub=='Reactive Astrocyte').sum()
            a=(sub=='Astrocyte').sum()
            denom=r+a
            out[f'f{um}']=float(r/denom) if denom>=10 else np.nan
            out[f'n{um}']=int(denom); out[f'r{um}']=int(r); out[f'a{um}']=int(a)
    rows.append(out)
df=pd.DataFrame(rows)
print(df[['donor_id','ca1_pyr','ca1_lineage','n20','f20','n30','f30']].to_string(index=False))
for f in ['f20','f30']:
    print('\n',f)
    print('n analyzable',df[f].notna().sum())
    print('partial',partial_correlation(df.rename(columns={'outcome':'slope_zmem0'}),feature_col=f,outcome_col='slope_zmem0',confounds=['max_age_vis','braak_numeric','cerad_ordinal','sex_binary']))
    print('loo',residualized_loo_predictive_correlation(df.rename(columns={'outcome':'slope_zmem0'}),feature_col=f,outcome_col='slope_zmem0',confounds=['max_age_vis','braak_numeric','cerad_ordinal','sex_binary']))
PY
