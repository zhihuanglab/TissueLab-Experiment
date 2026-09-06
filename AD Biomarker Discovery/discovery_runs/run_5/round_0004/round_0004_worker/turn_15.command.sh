set -euo pipefail
python - <<'PY'
import json, pandas as pd
res=json.load(open('/scratch/results.json'))
for r in res['ranked_variations']:
    print('\n',r['name'])
    for k in ['partial_r','p_value','ci_lo','ci_hi','selection_score','loo_predictive_r','gap','adjusted_score','bootstrap_median_partial_r','bootstrap_sign_consistency','loo_unstable_count','loo_max_shift']:
        print(k, r[k])
table=pd.read_csv('/scratch/donor_feature_table.csv')
print('\nFeature summary:')
for c in [col for col in table.columns if 'area_' in col and 'um2' in col and 'iqr' not in col and 'mean' not in col]:
    print(c, table[c].describe())
print('\ncounts summary', table['ca1_peripyramidal_reactive_astro_count_30um'].describe())
PY
