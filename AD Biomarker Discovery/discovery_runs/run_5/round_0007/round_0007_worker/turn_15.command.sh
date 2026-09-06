set -euo pipefail
python - <<'PY'
import json, pandas as pd, numpy as np
res=json.load(open('/scratch/results.json'))
df=pd.DataFrame(res['per_donor_loo'])
df['abs_error']=(df['outcome']-df['predicted']).abs()
print(df.sort_values('abs_error', ascending=False)[['donor_id','outcome','predicted',res['feature_column'],'n_ca1_lymphocytes','n_peripyramidal_reactive_astrocytes','abs_error']].head(10).to_string(index=False))
print('\nNonzero feature donors:')
print(df[df[res['feature_column']]>0][['donor_id','outcome',res['feature_column'],'n_ca1_lymphocytes']].to_string(index=False))
PY
