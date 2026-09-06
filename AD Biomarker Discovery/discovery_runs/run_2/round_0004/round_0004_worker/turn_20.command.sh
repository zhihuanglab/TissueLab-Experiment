set -euo pipefail
python - <<'PY'
import json, pandas as pd, numpy as np
res=json.loads(open('/scratch/results.json').read())
best=[r for r in res['ranked_variations'] if r['variation_name']==res['best_variation']][0]
df=pd.DataFrame(best['per_donor_loo'])
df['err']=df['actual_residualized_outcome']-df['predicted']
df['abs_err']=df['err'].abs()
top=df.sort_values('abs_err',ascending=False).head(10)
print(top[['donor_id','sex',res['feature_column'],'predicted','actual_residualized_outcome','err','ca1_astrocyte_count','ca1_reactive_astrocyte_count']].to_string(index=False))
print('\nTop10 sex counts:\n', top['sex'].value_counts().to_string())
print('\nPositive contrast donors:')
print(df[df[res['feature_column']]>0][['donor_id',res['feature_column'],'predicted','actual_residualized_outcome','sex','ca1_astrocyte_count']].to_string(index=False))
PY
