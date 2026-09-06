set -euo pipefail
python - <<'PY'
import json, pandas as pd
with open('/scratch/results.json') as f:
    obj=json.load(f)
feat=obj['feature_column']
df=pd.DataFrame(obj['winner_loo_table'])
df['error']=df['outcome']-df['predicted']
df['abs_error']=df['error'].abs()
print(df.sort_values('abs_error', ascending=False).head(10).to_string(index=False))
PY
