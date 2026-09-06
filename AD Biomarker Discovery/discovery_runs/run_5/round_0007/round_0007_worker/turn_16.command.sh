set -euo pipefail
python - <<'PY'
import pandas as pd
df=pd.read_csv('/scratch/donor_feature_table.csv')
for col in [c for c in df.columns if 'lymphocyte_contact' in c]:
    nz=(df[col]>0).sum()
    print(col, 'nonzero', nz, 'max', df[col].max(), 'mean', df[col].mean())
PY
