set -euo pipefail
python - <<'PY'
import json, pandas as pd, numpy as np
from pathlib import Path

results = json.loads(Path('/scratch/results.json').read_text())
df = pd.read_csv('/scratch/donor_feature_table.csv')
best_col = results['feature_column']
loo = pd.DataFrame(results['per_donor_loo'])
loo['abs_error'] = (loo['outcome'] - loo['predicted']).abs()
print('best_col', best_col)
print('nonzero count', int((df[best_col] > 0).sum()), 'of', len(df))
print('feature summary')
print(df[best_col].describe().to_string())
print('\nTop nonzero donors:')
print(df.loc[df[best_col] > 0, ['donor_id', best_col, 'n_ra_exposed_pyramidal_ca1', 'n_lymphocyte_ca1']].sort_values(best_col, ascending=False).to_string(index=False))
print('\nLargest absolute LOO errors:')
print(loo.sort_values('abs_error', ascending=False).head(10).to_string(index=False))
PY
