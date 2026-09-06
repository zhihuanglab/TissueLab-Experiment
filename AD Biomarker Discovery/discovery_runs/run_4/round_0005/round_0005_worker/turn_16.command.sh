set -euo pipefail
python - <<'PY'
import json, pandas as pd, sys
sys.path.insert(0,'/shared/lib')
from shared_analysis.sea_ad_lfb import load_training_cohort
res=json.load(open('/scratch/results.json'))
loo=pd.DataFrame(res['loo_rows'])
loo['abs_err']=(loo['predicted']-loo['outcome']).abs()
cohort=load_training_cohort('/data')
df=pd.read_csv('/scratch/donor_feature_table.csv')
winner=res['feature_column']
merged=cohort.merge(df[['donor_id',winner,'ca1_reactive_n','ca1_astro_lineage_n','reactive_purity_r60um_valid_centers']], on='donor_id').merge(loo[['donor_id','abs_err','predicted']], on='donor_id')
print(merged.sort_values('abs_err', ascending=False).head(10)[['donor_id','slope_zmem0','predicted','abs_err','sex','cognitive_status','braak_numeric','cerad_ordinal',winner,'ca1_reactive_n','reactive_purity_r60um_valid_centers']].to_string(index=False))
PY
