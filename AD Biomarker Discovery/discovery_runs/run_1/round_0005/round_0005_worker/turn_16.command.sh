set -e
echo '--- results.json ---'
python - <<'PY'
import json
with open('/scratch/results.json') as f:
    obj=json.load(f)
for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','gap_penalty','adjusted_score']:
    print(k, ':', obj.get(k))
print('ranked n', len(obj.get('ranked_variations', [])))
print('first ranked keys', obj['ranked_variations'][0].keys())
PY
echo '--- report.md ---'
sed -n '1,260p' /scratch/report.md
echo '--- donor_feature_table columns ---'
python - <<'PY'
import pandas as pd
df=pd.read_csv('/scratch/donor_feature_table.csv')
print(df.columns.tolist())
print(df.head(2).to_string(index=False))
PY
