set -euo pipefail
python - <<'PY'
import json
import pandas as pd
from pathlib import Path

results=json.loads(Path('/scratch/results.json').read_text())
best=results['best_variation']
best_row=[r for r in results['ranked_variations'] if r['name']==best][0]
count_col=best_row['selected_count_column']
df=pd.read_csv('/scratch/donor_feature_table.csv')
print('feature summary')
print(df[[results['feature_column'], count_col, 'ca1_pyramidal_neuron_count', 'ca1_reactive_astrocyte_count']].describe().to_string())

lines=Path('/scratch/run_stdout.txt').read_text().splitlines()
start=lines.index('PER-DONOR (LOO):')
rows=[]
for line in lines[start+2:]:
    if not line.strip():
        continue
    parts=line.split()
    rows.append({
        'donor_id': parts[0],
        'outcome': float(parts[1]),
        'predicted': float(parts[2]),
        'feature': float(parts[3]),
    })
pred=pd.DataFrame(rows)
pred['abs_error']=(pred['outcome']-pred['predicted']).abs()
cohort=pd.read_csv('/data/training_cohort.csv')[['donor_id','braak_numeric','cerad_ordinal','sex','max_age_vis']]
pred=pred.merge(cohort,on='donor_id',how='left').sort_values('abs_error', ascending=False)
print('\nTop errors')
print(pred.head(10).to_string(index=False))
print('\nError by sex')
print(pred.groupby('sex')['abs_error'].agg(['mean','median','count']).to_string())
print('\nError corr with pathology/age')
print(pred[['abs_error','braak_numeric','cerad_ordinal','max_age_vis','feature']].corr(numeric_only=True).loc['abs_error'].to_string())
PY
