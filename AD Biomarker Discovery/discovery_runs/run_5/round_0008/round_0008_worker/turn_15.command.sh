set -euo pipefail
python - <<'PY'
import pandas as pd
df=pd.read_csv('/scratch/donor_feature_table.csv')
col='ca1_pyramidal_reactive_astro_lymphocyte_triad_fraction_25um'
print(df[['donor_id','slope_zmem0',col,'n_ca1_pyramidal_neurons','n_ca1_reactive_astrocytes','n_ca1_lymphocytes']].sort_values(col, ascending=False).head(12).to_string(index=False))
print('\nWorst outcomes:')
print(df[['donor_id','slope_zmem0',col,'n_ca1_lymphocytes']].sort_values('slope_zmem0').head(12).to_string(index=False))
PY
