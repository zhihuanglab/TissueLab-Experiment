set -euo pipefail
python - <<'PY'
import pandas as pd
df=pd.read_csv('/scratch/donor_feature_table.csv')
print(df.columns.tolist())
print(df.head(2).to_string(index=False))
PY
