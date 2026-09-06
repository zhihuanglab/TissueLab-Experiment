set -euo pipefail
python - <<'PY'
import json
import pandas as pd
from pathlib import Path

results=json.loads(Path('/scratch/results.json').read_text())
df=pd.read_csv('/scratch/donor_feature_table.csv')
best=results['best_variation']
feature=results['feature_column']
best_row=[r for r in results['ranked_variations'] if r['name']==best][0]
print('best', best_row)

# counts summary
count_col=best_row['selected_count_column']
print('\ncount summary')
print(df[count_col].describe().to_string())

# merge run stdout loo predictions
lines=Path('/scratch/run_stdout.txt').read_text().splitlines()
start=lines.index('PER-DONOR (LOO):')
header=lines[start+1]
rows=[]
for line in lines[start+2:]:
    if not line.strip(): continue
    parts=line.split()
    donor=parts[0]
    outcome=float(parts[1]); pred=float(parts[2]); feat=float(parts[3])
    rows.append({'donor_id':donor,'outcome':outcome,'predicted':pred,'feature':feat,'abs_error':abs(outcome-pred)})
pred=pd.DataFrame(rows)
cohort=pd.read_csv('/data/training_cohort.csv')[['donor_id','cognitive_status','overall_ad_neuropath_change','braak_numeric','cerad_ordinal','sex']]
pred=pred.merge(cohort,on='donor_id',how='left').sort_values('abs_error', ascending=False)
print('\nTop abs errors:')
print(pred.head(10).to_string(index=False))

print('\nMean abs error by cognitive_status')
print(pred.groupby('cognitive_status')['abs_error'].mean().to_string())
print('\nMean abs error by AD change')
print(pred.groupby('overall_ad_neuropath_change')['abs_error'].mean().to_string())
PY
