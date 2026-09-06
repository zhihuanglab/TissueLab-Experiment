set -euo pipefail
python - <<'PY'
import json, pandas as pd, numpy as np
payload=json.load(open('/scratch/results.json'))
print(payload.keys())
print('feature_column', payload['feature_column'])
df=pd.DataFrame(payload['per_donor'])
print(df[['donor_id','outcome','predicted','abs_error',payload['feature_column'],payload['numerator_column'] if 'numerator_column' in payload else df.columns[0]] if False else ['donor_id']].head())
print('top abs errors')
print(df[['donor_id','outcome','predicted','abs_error',payload['feature_column']]].sort_values('abs_error', ascending=False).head(10).to_string(index=False))
PY
