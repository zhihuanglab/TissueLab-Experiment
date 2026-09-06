set -euo pipefail
python - <<'PY'
import sys, json, math
sys.path.append('/shared/lib')
from shared_analysis import load_training_cohort, build_cell_table
from pathlib import Path
import pandas as pd
data_root=Path('/data')
cohort=load_training_cohort(data_root)
rows=[]
for _,r in cohort.iterrows():
    slide=data_root / r['slide_name']
    df=build_cell_table(slide, include_regions=True, include_geometry=False)
    sub=df[df['region']=='CA1']
    counts=sub['cell_type'].value_counts()
    rows.append({
        'donor_id':r['donor_id'],
        'n_ca1_pyr': int(counts.get('Pyramidal Neuron',0)),
        'n_ca1_ra': int(counts.get('Reactive Astrocyte',0)),
        'n_ca1_lym': int(counts.get('Lymphocyte',0)),
        'n_ca1_total': int(len(sub)),
    })
out=pd.DataFrame(rows)
print(out.describe(include='all').to_string())
print('\nMin counts donors:')
print(out.sort_values(['n_ca1_pyr','n_ca1_ra','n_ca1_lym']).head(10).to_string(index=False))
print('\nZero-count summary:')
print((out[['n_ca1_pyr','n_ca1_ra','n_ca1_lym']]==0).sum().to_string())
PY
