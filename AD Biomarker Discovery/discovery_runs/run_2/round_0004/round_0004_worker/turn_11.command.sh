set -euo pipefail
python - <<'PY'
import sys, warnings, numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0,'/shared/lib')
from shared_analysis import load_training_cohort, build_cell_table
warnings.filterwarnings('ignore')
cohort = load_training_cohort('/data')
rows=[]
for _, row in cohort.iterrows():
    cells=build_cell_table(Path('/data')/row.slide_name, include_regions=True, include_geometry=False)
    ca1=cells[cells['region']=='CA1']
    vc=ca1['cell_type'].value_counts()
    rows.append({
        'donor_id':row.donor_id,
        'reactive':int(vc.get('Reactive Astrocyte',0)),
        'astro':int(vc.get('Astrocyte',0)),
        'ca1_total':int(len(ca1)),
    })
df=pd.DataFrame(rows)
print(df[['reactive','astro']].describe().to_string())
print('donors with either class <20:', df[(df.reactive<20)|(df.astro<20)][['donor_id','reactive','astro']].to_string(index=False))
PY
