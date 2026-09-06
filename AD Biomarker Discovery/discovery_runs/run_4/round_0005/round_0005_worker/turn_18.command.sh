set -euo pipefail
python - <<'PY'
import json, pandas as pd
res=json.load(open('/scratch/results.json'))
loo=pd.DataFrame(res['loo_rows'])
loo['abs_err']=(loo['predicted']-loo['outcome']).abs()
cohort=pd.read_csv('/data/training_cohort.csv')
df=pd.read_csv('/scratch/donor_feature_table.csv')
winner=res['feature_column']
merged=cohort.merge(df[['donor_id',winner,'ca1_reactive_n','ca1_astro_lineage_n','reactive_purity_r60um_valid_centers','reactive_purity_r40um_valid_centers']], on='donor_id').merge(loo[['donor_id','abs_err','predicted']], on='donor_id')
print(merged.sort_values('abs_err', ascending=False).head(10)[['donor_id','slope_zmem0','predicted','abs_err','sex','braak_numeric','cerad_ordinal',winner,'ca1_reactive_n','reactive_purity_r60um_valid_centers','reactive_purity_r40um_valid_centers']].to_string(index=False))
PY
