set -euo pipefail
python - <<'PY'
import json, pandas as pd
res=json.load(open('/scratch/results.json'))
print('best', res['best_variation'], 'feature_column', res['feature_column'])
print('top metrics:', {k: res[k] for k in ['n_total','n_analyzable','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score','p_value','ci_lo','ci_hi','bootstrap_median_partial_r','bootstrap_sign_consistency','loo_unstable_count','loo_max_shift']})
best=res['ranked_variations'][0]
df=pd.DataFrame(best['per_donor_loo'])
cohort=pd.read_csv('/data/training_cohort.csv')
merged=df.merge(cohort[['donor_id','sex','cognitive_status','braak_numeric','cerad_ordinal','overall_ad_neuropath_change']], on='donor_id', how='left')
print('\nTop abs errors:')
print(merged.sort_values('abs_error', ascending=False).head(10).to_string(index=False))
print('\nVariation ranking extra:')
for r in res['ranked_variations']:
    print(r['name'], r['feature_column'], r['partial_r'], r['loo_predictive_r'], r['p_value'], r['bootstrap_median_partial_r'], r['bootstrap_sign_consistency'])
PY
