set -euo pipefail
python - <<'PY'
import pandas as pd
df=pd.read_csv('/scratch/donor_feature_table.csv')
feat='ca1_reactive_niche_pyramidal_corpora_any_fraction_r96px'
print('nonzero donors:', df.loc[df[feat]>0, ['donor_id', feat, 'ca1_pyramidal_n','ca1_reactive_astrocyte_n','ca1_corpora_amylacea_n','reactive_niche_pyramidal_n']].to_string(index=False))
print('summary', df[feat].describe())
PY
