set -e
python - <<'PY'
import pandas as pd
df=pd.read_csv('/scratch/donor_feature_table.csv')
for c in [col for col in df.columns if 'neighbor_fraction' in col]:
    nz=(df[c]>0).sum()
    print(c, 'nonzero', nz, 'max', df[c].max(), 'median', df[c].median())
PY
