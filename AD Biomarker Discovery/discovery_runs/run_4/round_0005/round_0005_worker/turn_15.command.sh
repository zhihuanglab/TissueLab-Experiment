set -euo pipefail
python - <<'PY'
import json, pandas as pd, numpy as np
res=json.load(open('/scratch/results.json'))
loo=pd.DataFrame(res['loo_rows'])
loo['abs_err']=(loo['predicted']-loo['outcome']).abs()
print('Worst abs errors:')
print(loo.sort_values('abs_err', ascending=False).head(8)[['donor_id','outcome','predicted','abs_err']].to_string(index=False))
df=pd.read_csv('/scratch/donor_feature_table.csv')
winner=res['feature_column']
print('\nWinner feature summary:')
print(df[winner].describe().to_string())
print('\nCounts summary:')
for c in ['ca1_astro_lineage_n','ca1_reactive_n','reactive_purity_r40um_valid_centers','reactive_purity_r60um_valid_centers']:
    print(c, df[c].describe().to_string())
print('\nWorst donor covariates:')
print(df.merge(loo[['donor_id','abs_err']], on='donor_id').sort_values('abs_err', ascending=False).head(8)[['donor_id','slope_zmem0','abs_err','sex','cognitive_status','braak_numeric','cerad_ordinal',winner,'ca1_reactive_n']].to_string(index=False))
PY
