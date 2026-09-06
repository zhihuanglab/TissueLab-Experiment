set -euo pipefail
python - <<'PY'
import pandas as pd
df=pd.read_csv('/scratch/donor_feature_table.csv')
print('rows', len(df))
print('missing slide_mpp', df['slide_mpp'].isna().sum())
print(df['slide_mpp'].describe())
PY
