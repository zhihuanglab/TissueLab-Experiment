set -euo pipefail
python - <<'PY'
import json, pandas as pd
with open('/scratch/results.json') as f:
    res=json.load(f)
bestcol=res['feature_column']
per=pd.DataFrame(res.get('ranked_variations',[]))
print('bestcol', bestcol)
# parse per-donor loo from script? not in json. recompute from stdout unavailable.
# use donor_feature_table + manual hardcoded from previous stdout? easier parse maybe from report? instead inspect donor_feature_table only.
df=pd.read_csv('/scratch/donor_feature_table.csv')
print(df.columns.tolist())
print(df[['donor_id',bestcol,'ca1_pyramidal_n','ca1_reactive_astrocyte_n','ca1_oligodendrocyte_n','sex','braak_numeric','cerad_ordinal']].head())
PY
