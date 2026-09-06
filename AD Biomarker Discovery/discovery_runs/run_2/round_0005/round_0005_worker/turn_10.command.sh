set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0,'/shared/lib')
from pathlib import Path
from shared_analysis import load_training_cohort, build_cell_table
import pandas as pd, numpy as np

cohort = load_training_cohort('/data')
rows=[]
for _, r in cohort.iterrows():
    slide = Path('/data')/r['slide_name']
    cells = build_cell_table(slide, include_regions=True, include_geometry=False)
    ca1 = cells[cells['region']=='CA1']
    rows.append({
        'donor_id': r['donor_id'],
        'n_ca1': len(ca1),
        'n_reactive': int((ca1['cell_type']=='Reactive Astrocyte').sum()),
        'n_astro': int((ca1['cell_type']=='Astrocyte').sum()),
    })
df = pd.DataFrame(rows)
print(df.describe().T[['min','25%','50%','75%','max']])
print('\nlow reactive:')
print(df.nsmallest(10,'n_reactive')[['donor_id','n_reactive','n_astro','n_ca1']].to_string(index=False))
print('\nlow astro:')
print(df.nsmallest(10,'n_astro')[['donor_id','n_reactive','n_astro','n_ca1']].to_string(index=False))
PY
