set -euo pipefail
python - <<'PY'
import json, pandas as pd, numpy as np
with open('/scratch/results.json') as f:
    res=json.load(f)
print('best', res['best_variation'], res['feature_column'])
ranked = {r['variation_name']: r for r in res['ranked_variations']}
for name,row in ranked.items():
    per = pd.DataFrame(row['per_donor_loo'])
    print('\nVAR', name)
    print('top abs errors:')
    print(per[['donor_id','outcome','predicted','error','abs_error',row['feature_column']]].sort_values('abs_error', ascending=False).head(6).to_string(index=False))
    if 'unstable_donors' in row:
        print('unstable donors:', [d['donor_id'] for d in row['unstable_donors']])
feat = pd.read_csv('/scratch/donor_feature_table.csv')
best_col=res['feature_column']
cols=['donor_id','slope_zmem0','sex','braak_numeric','cerad_ordinal','cognitive_status','overall_ad_neuropath_change',best_col,'ca1_classified','ca2_classified']
print('\nTop high feature')
print(feat[cols].sort_values(best_col, ascending=False).head(8).to_string(index=False))
print('\nTop low feature')
print(feat[cols].sort_values(best_col, ascending=True).head(8).to_string(index=False))
print('\nOutcome summary', feat['slope_zmem0'].describe().to_string())
print('Feature summary', feat[best_col].describe().to_string())
print('\nCorr raw:', feat[[best_col,'slope_zmem0']].corr().iloc[0,1])
PY
