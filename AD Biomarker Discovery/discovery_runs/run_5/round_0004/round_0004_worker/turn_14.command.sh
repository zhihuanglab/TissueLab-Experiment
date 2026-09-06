set -euo pipefail
python - <<'PY'
import pandas as pd, json
cohort=pd.read_csv('/data/training_cohort.csv')
print(cohort.columns.tolist())
res=json.load(open('/scratch/results.json'))
best=res['ranked_variations'][0]
df=pd.DataFrame(best['per_donor_loo'])
avail=[c for c in ['donor_id','sex','braak_numeric','cerad_ordinal','braak_label','cerad_label','overall_ad_neuropath_change','cogdx','ADNC'] if c in cohort.columns]
merged=df.merge(cohort[avail], on='donor_id', how='left')
print('\nTop abs errors:')
print(merged.sort_values('abs_error', ascending=False).head(10).to_string(index=False))
print('\nUnstable donors', best['unstable_donors'])
PY
