set -euo pipefail
python - <<'PY'
import json, pandas as pd, numpy as np
from pathlib import Path
res=json.loads(Path('/scratch/results.json').read_text())
print('best', res['best_variation'], res['feature_column'])
print('ranked', [(r['variation_name'], r['partial_r'], r['loo_predictive_r']) for r in res['ranked_variations']])
txt=Path('/scratch/result_stdout.txt').read_text().splitlines()
rows=[]
start=False
for line in txt:
    if line.startswith('  donor_id  outcome'):
        start=True
        continue
    if start and line.strip():
        parts=line.split()
        # donor_id outcome predicted feature pyr react lymph
        rows.append({
            'donor_id':parts[0],
            'outcome':float(parts[1]),
            'predicted':float(parts[2]),
            'feature':float(parts[3]),
            'pyr':int(parts[4]),
            'react':int(parts[5]),
            'lymph':int(parts[6]),
        })
df=pd.DataFrame(rows)
df['error']=df['predicted']-df['outcome']
df['abs_error']=df['error'].abs()
print('\nTop abs errors:')
print(df.sort_values('abs_error', ascending=False).head(10).to_string(index=False))
print('\nZero-feature donors:', (df['feature']==0).sum(), 'of', len(df))
print(df.groupby(df['feature']==0)['outcome'].agg(['count','mean']))
print('\nHigh-feature donors:')
print(df.sort_values('feature', ascending=False).head(8).to_string(index=False))
PY
