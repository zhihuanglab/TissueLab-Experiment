set -euo pipefail
python - <<'PY'
import json, pandas as pd, numpy as np
from pathlib import Path
res=json.loads(Path('/scratch/results.json').read_text())
best=[r for r in res['ranked_variations'] if r['variation_name']==res['best_variation']][0]
df=pd.DataFrame(best['per_donor_loo'])
df['abs_err']=(df['actual_residualized_outcome']-df['predicted']).abs()
print('Top abs errors:')
print(df.sort_values('abs_err', ascending=False)[['donor_id','outcome','predicted','actual_residualized_outcome',res['feature_column'],'ca1_reactive_astrocyte_count','ca1_astrocyte_count','sex']].head(10).to_string(index=False))
print('\nPrediction sign agreement rate', np.mean(np.sign(df['predicted'])==np.sign(df['actual_residualized_outcome'])))
print('Predicted summary:', df['predicted'].describe().to_string())
print('Feature summary:', df[res['feature_column']].describe().to_string())
print('\nCorrelations:')
for col in [res['feature_column'],'reactive_median_log1p_area_ca1','astrocyte_median_log1p_area_ca1','ca1_reactive_astrocyte_count','ca1_astrocyte_count']:
    print(col, df[col].corr(df['outcome']))
PY
