set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0, '/shared/lib')
from shared_analysis import load_training_cohort, build_cell_table
import pandas as pd, numpy as np
cohort = load_training_cohort('/data')
print('cohort rows', len(cohort))
slide = '/data/' + cohort.iloc[0]['slide_name']
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print('slide', slide)
print('cell types top 15:')
print(cells['cell_type'].value_counts().head(15).to_string())
print('regions counts:')
print(cells['region'].value_counts(dropna=False).to_string())
for reg in ['CA1','CA2']:
    sub = cells[cells['region']==reg]
    classified = sub[~sub['cell_type'].str.lower().isin(['unknown','unlabeled'])]
    print(reg, 'n', len(sub), 'classified', len(classified), 'pyr_frac', (classified['cell_type']=='Pyramidal Neuron').mean())
PY
