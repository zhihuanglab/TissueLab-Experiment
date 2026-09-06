set -euo pipefail
python - <<'PY'
import json, pandas as pd
obj=json.load(open('/scratch/results.json'))
df=pd.DataFrame(obj['per_donor_loo'])
df['abs_error']=(df['predicted']-df['slope_zmem0']).abs()
print(df.sort_values('abs_error', ascending=False).head(10).to_string(index=False))
PY
