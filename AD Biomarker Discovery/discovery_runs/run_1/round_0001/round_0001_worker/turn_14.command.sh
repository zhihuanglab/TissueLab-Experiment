set -euo pipefail
python - <<'PY'
import pandas as pd
feat = pd.read_csv('/scratch/donor_feature_table.csv')
print(feat.columns.tolist())
print(feat[['donor_id','slope_zmem0','sex','braak_numeric','cerad_ordinal','ca1_pyramidal_fraction','ca1_classified','ca2_classified']].sort_values('ca1_pyramidal_fraction', ascending=False).head(8).to_string(index=False))
print('\nLowest')
print(feat[['donor_id','slope_zmem0','sex','braak_numeric','cerad_ordinal','ca1_pyramidal_fraction','ca1_classified','ca2_classified']].sort_values('ca1_pyramidal_fraction', ascending=True).head(8).to_string(index=False))
PY
