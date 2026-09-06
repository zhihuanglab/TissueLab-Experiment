set -euo pipefail
python - <<'PY'
import pandas as pd, json
feat = pd.read_csv('/scratch/donor_feature_table.csv')
with open('/scratch/results.json') as f: res = json.load(f)
best = [r for r in res['ranked_variations'] if r['variation_name']==res['best_variation']][0]
per = pd.DataFrame(best['per_donor_loo'])
meta = feat[['donor_id','sex','braak_numeric','cerad_ordinal','max_age_vis','ca1_pyramidal_fraction','ca1_ca2_pyramidal_fraction','ca1_classified','ca2_classified']].merge(per[['donor_id','outcome','predicted','error','abs_error']], on='donor_id')
print(meta.sort_values('abs_error', ascending=False).head(10).to_string(index=False))
PY
